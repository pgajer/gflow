#' Create Uniform 2D Grid
#'
#' @description Creates a uniform 2D grid over \code{[0,1]^2} for function evaluation.
#'
#' @param axis.n.pts Integer number of points along each axis
#' @return List containing:
#'   \item{x}{Numeric vector of x coordinates}
#'   \item{y}{Numeric vector of y coordinates}
#'   \item{grid}{Data frame with all (x,y) coordinate pairs}
#'
#' @examples
#' grid <- create.grid(50)
#' str(grid)
#' @export
create.grid <- function(axis.n.pts) {
    ## Creates uniform 2D grid over [0,1]^2
    ## Returns: list containing x, y vectors and grid data.frame
    x <- seq(0, 1, length.out = axis.n.pts)
    y <- seq(0, 1, length.out = axis.n.pts)
    grid <- expand.grid(x = x, y = y)
    return(list(x = x, y = y, grid = grid))
}

#' Evaluate Function on Grid as Vector
#'
#' @description Evaluates a bivariate function at each point of a grid, returning results as a vector.
#'
#' @param f Function taking two arguments (x, y)
#' @param grid List containing grid data (as created by create.grid)
#' @return List with components:
#'   \item{y.smooth}{Numeric vector of function values}
#'   \item{X}{Matrix of grid coordinates}
#'
#' @examples
#' grid <- create.grid(50)
#' f <- function(x, y) x^2 + y^2
#' values <- evaluate.function.on.grid.as.vector(f, grid)
#' @export
evaluate.function.on.grid.as.vector <- function(f, grid) {
    # Convert grid to matrix of coordinates
    X <- as.matrix(grid$grid)

    # Initialize output vector with length equal to number of rows in X
    y.smooth <- numeric(nrow(X))

    # Evaluate function at each point in X
    for(i in 1:nrow(X)) {
        y.smooth[i] <- f(X[i, 1], X[i, 2])
    }

    list(y.smooth = y.smooth,
         X = X)
}

#' Evaluate Function on Grid
#'
#' @description Evaluates a bivariate function on a uniform grid.
#'
#' @param f Function taking two arguments (x, y)
#' @param grid List containing grid data (as created by create.grid)
#' @return Matrix of function values
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(50)
#' f <- function(x, y) sin(x) * cos(y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' image(grid$x, grid$y, f.grid)
#' }
#' @export
evaluate.function.on.grid <- function(f, grid) {
    ## Evaluates scalar function f on grid
    ## Returns: matrix of function values
    n <- length(grid$x)
    f.grid <- matrix(0, n, n)

    for (i in 1:n) {
        for (j in 1:n) {
            f.grid[i,j] <- f(grid$x[i], grid$y[j])
        }
    }

    return(f.grid)
}

## Neighbor computation functions
get.neighbors <- function(i, j, n) {
    ## Get valid grid neighbors for point (i,j) in nxn grid
    ## Returns: matrix of neighbor coordinates
    neighbors <- matrix(c(
        i+1, j,    # right
        i-1, j,    # left
        i, j+1,    # up
        i, j-1     # down
    ), ncol = 2, byrow = TRUE)

    ## Filter out invalid neighbors
    valid <- neighbors[,1] >= 1 & neighbors[,1] <= n &
             neighbors[,2] >= 1 & neighbors[,2] <= n

    return(neighbors[valid, , drop = FALSE])
}

## Critical point detection
is.local.maximum <- function(i, j, f.grid) {
    ## Check if (i,j) is a local maximum
    n <- nrow(f.grid)
    neighbors <- get.neighbors(i, j, n)

    if (nrow(neighbors) == 0) return(TRUE)  # Edge case

    val <- f.grid[i,j]
    for (k in 1:nrow(neighbors)) {
        if (f.grid[neighbors[k,1], neighbors[k,2]] > val) {
            return(FALSE)
        }
    }
    return(TRUE)
}

is.local.minimum <- function(i, j, f.grid) {
    ## Check if (i,j) is a local minimum
    n <- nrow(f.grid)
    neighbors <- get.neighbors(i, j, n)

    if (nrow(neighbors) == 0) return(TRUE)

    val <- f.grid[i,j]
    for (k in 1:nrow(neighbors)) {
        if (f.grid[neighbors[k,1], neighbors[k,2]] < val) {
            return(FALSE)
        }
    }
    return(TRUE)
}

#' Find Critical Points in Function Grid
#'
#' @description Identifies local maxima and minima in a gridded function.
#'
#' @param f.grid Matrix of function values
#' @return List containing:
#'   \item{maxima}{Matrix of (i,j) coordinates for local maxima}
#'   \item{minima}{Matrix of (i,j) coordinates for local minima}
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(50)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' critical <- find.critical.points(f.grid)
#' print(critical)
#' }
#' @export
find.critical.points <- function(f.grid) {
    ## Find all critical points in the grid
    n <- nrow(f.grid)
    maxima <- matrix(nrow = 0, ncol = 2)
    minima <- matrix(nrow = 0, ncol = 2)

    for (i in 1:n) {
        for (j in 1:n) {
            if (is.local.maximum(i, j, f.grid)) {
                maxima <- rbind(maxima, c(i, j))
            }
            if (is.local.minimum(i, j, f.grid)) {
                minima <- rbind(minima, c(i, j))
            }
        }
    }

    return(list(maxima = maxima, minima = minima))
}

## Gradient flow computation
find.steepest.ascent.neighbor <- function(i, j, f.grid) {
    ## Find neighbor in direction of steepest ascent
    n <- nrow(f.grid)
    neighbors <- get.neighbors(i, j, n)

    if (nrow(neighbors) == 0) return(c(i, j))

    max.val <- f.grid[i, j]
    max.idx <- c(i, j)

    for (k in 1:nrow(neighbors)) {
        ni <- neighbors[k, 1]
        nj <- neighbors[k, 2]
        if (f.grid[ni, nj] > max.val) {
            max.val <- f.grid[ni, nj]
            max.idx <- c(ni, nj)
        }
    }

    return(max.idx)
}

find.steepest.descent.neighbor <- function(i, j, f.grid) {
    ## Find neighbor in direction of steepest descent
    n <- nrow(f.grid)
    neighbors <- get.neighbors(i, j, n)

    if (nrow(neighbors) == 0) return(c(i, j))

    min.val <- f.grid[i, j]
    min.idx <- c(i, j)

    for (k in 1:nrow(neighbors)) {
        ni <- neighbors[k, 1]
        nj <- neighbors[k, 2]
        if (f.grid[ni, nj] < min.val) {
            min.val <- f.grid[ni, nj]
            min.idx <- c(ni, nj)
        }
    }

    return(min.idx)
}

#' Compute Gradient Trajectory
#'
#' @description Computes both ascending and descending gradient trajectories from a grid point.
#'
#' @param i Row index of starting point
#' @param j Column index of starting point
#' @param f.grid Matrix of function values
#' @return List containing:
#'   \item{ascending.path}{Matrix of (i,j) coordinates along ascending path}
#'   \item{descending.path}{Matrix of (i,j) coordinates along descending path}
#'   \item{local.max}{Coordinates of destination local maximum}
#'   \item{local.min}{Coordinates of destination local minimum}
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(30)
#' f <- function(x, y) -((x-0.5)^2 + (y-0.5)^2)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' traj <- compute.gradient.trajectory(15, 15, f.grid)
#' str(traj)
#' }
#' @export
compute.gradient.trajectory <- function(i, j, f.grid) {
    ## Compute full gradient trajectory from point (i,j)
    n <- nrow(f.grid)
    max.steps <- n * n  # Safety limit

    ## Ascending trajectory
    asc.path <- matrix(c(i, j), nrow = 1)
    curr <- c(i, j)

    for (step in 1:max.steps) {
        next.pos <- find.steepest.ascent.neighbor(curr[1], curr[2], f.grid)
        if (all(next.pos == curr)) break  # At local maximum

        curr <- next.pos
        asc.path <- rbind(asc.path, curr)
    }

    ## Descending trajectory
    desc.path <- matrix(c(i, j), nrow = 1)
    curr <- c(i, j)

    for (step in 1:max.steps) {
        next.pos <- find.steepest.descent.neighbor(curr[1], curr[2], f.grid)
        if (all(next.pos == curr)) break  # At local minimum

        curr <- next.pos
        desc.path <- rbind(desc.path, curr)
    }

    return(list(
        ascending.path = asc.path,
        descending.path = desc.path,
        local.max = asc.path[nrow(asc.path), ],
        local.min = desc.path[nrow(desc.path), ]
    ))
}

#' Compute Morse-Smale Cells
#'
#' @description Assigns each grid point to a Morse-Smale cell based on gradient flow.
#'
#' @param f.grid Matrix of function values
#' @return Matrix of cell labels
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(30)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' cells <- compute.morse.smale.cells(f.grid)
#' image(cells)
#' }
#' @export
compute.morse.smale.cells <- function(f.grid) {
    n <- nrow(f.grid)
    cell.labels <- matrix(0, n, n)

    ## Create unique label for each min-max pair
    label.counter <- 1
    label.map <- list()

    ## Compute gradient trajectory for each point and assign cell labels
    for (i in 1:n) {
        for (j in 1:n) {
            traj <- compute.gradient.trajectory(i, j, f.grid)

            ## Create key for min-max pair
            key <- paste(paste(traj$local.min, collapse=","),
                        paste(traj$local.max, collapse=","),
                        sep="|")

            ## Assign or create new label
            if (is.null(label.map[[key]])) {
                label.map[[key]] <- label.counter
                label.counter <- label.counter + 1
            }

            cell.labels[i,j] <- label.map[[key]]
        }
    }

    return(cell.labels)
}

#' Find Critical Points with Continuous Coordinates
#'
#' @description Identifies critical points (maxima, minima, saddles) using continuous gradient estimation.
#'
#' @param gradient Function returning gradient vector at (x,y)
#' @param grid Grid object created by create.grid
#' @return List containing matrices of critical point coordinates
#' @export
find.critical.points.continuous <- function(gradient, grid) {
    maxima <- matrix(nrow = 0, ncol = 2)
    minima <- matrix(nrow = 0, ncol = 2)
    saddles <- matrix(nrow = 0, ncol = 2)

    # Initialize eigenvalue storage
    attr(saddles, "eigenvalues") <- matrix(nrow = 0, ncol = 2)

    n <- length(grid$x)
    h <- 1e-5  # Step size for numerical Hessian

    for (i in 1:n) {
        for (j in 1:n) {
            x <- grid$x[i]
            y <- grid$y[j]

            # Check if gradient is near zero
            grad <- gradient(x, y)
            if (sqrt(sum(grad^2)) > 1e-3) next

            # Compute Hessian numerically
            H <- matrix(0, 2, 2)
            for (ii in 1:2) {
                for (jj in 1:2) {
                    H[ii, jj] <- numerical.second.derivative(gradient, x, y, ii, jj, h)
                }
            }

            # Force symmetry
            H <- 0.5 * (H + t(H))

            # Eigenvalue analysis
            eig <- eigen(H)
            lambda1 <- eig$values[1]
            lambda2 <- eig$values[2]

            # Classify critical point based on eigenvalues
            if (lambda1 > 0 && lambda2 > 0) {
                # Local minimum (positive definite)
                minima <- rbind(minima, c(x, y))
            } else if (lambda1 < 0 && lambda2 < 0) {
                # Local maximum (negative definite)
                maxima <- rbind(maxima, c(x, y))
            } else if (lambda1 * lambda2 < 0) {
                # Saddle point (eigenvalues have opposite signs)
                saddles <- rbind(saddles, c(x, y))

                # Store eigenvalues with each point for additional analysis
                current_eigenvalues <- attr(saddles, "eigenvalues")
                if (is.null(current_eigenvalues)) {
                    attr(saddles, "eigenvalues") <- matrix(c(lambda1, lambda2), nrow = 1)
                } else {
                    attr(saddles, "eigenvalues") <- rbind(current_eigenvalues, c(lambda1, lambda2))
                }
            }
        }
    }

    # Return all types of critical points
    return(list(
        maxima = maxima,
        minima = minima,
        saddles = saddles
    ))
}

# Helper function to compute numerical second derivatives
numerical.second.derivative <- function(gradient, x, y, i, j, h) {
    if (i == 1 && j == 1) {
        # ∂²f/∂x²
        g1 <- gradient(x + h, y)[1]
        g2 <- gradient(x - h, y)[1]
        return((g1 - g2) / (2 * h))
    } else if (i == 2 && j == 2) {
        # ∂²f/∂y²
        g1 <- gradient(x, y + h)[2]
        g2 <- gradient(x, y - h)[2]
        return((g1 - g2) / (2 * h))
    } else {
        # ∂²f/∂x∂y = ∂²f/∂y∂x
        g1 <- gradient(x + h, y + h)[1]
        g2 <- gradient(x - h, y - h)[1]
        return((g1 - g2) / (4 * h))
    }
}

compute.morse.smale.cells1 <- function(f.grid) {
    n <- nrow(f.grid)
    cell.labels <- matrix(0, n, n)
    critical.points <- find.critical.points(f.grid)

    ## Create unique label for each min-max pair
    label.counter <- 1
    label.map <- list()

    ## Compute gradient trajectory for each point and assign cell labels
    for (i in 1:n) {
        for (j in 1:n) {
            traj <- compute.gradient.trajectory(i, j, f.grid)

            ## Create key for min-max pair
            key <- paste(paste(traj$local.min, collapse=","),
                        paste(traj$local.max, collapse=","),
                        sep="|")

            ## Assign or create new label
            if (is.null(label.map[[key]])) {
                label.map[[key]] <- label.counter
                label.counter <- label.counter + 1
            }

            cell.labels[i,j] <- label.map[[key]]
        }
    }

    return(cell.labels)
}

## Visualization functions

#' Plot Gradient Trajectories
#'
#' @description Visualizes gradient flow trajectories from sample points on a contour plot.
#'
#' @param grid Grid object from create.grid
#' @param f.grid Matrix of function values
#' @param sample.points Matrix of (i,j) indices for sample trajectories (NULL for automatic)
#' @param spacing Integer spacing between sample points (default 10)
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(50)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' gradient.trajectories.plot(grid, f.grid, spacing = 15)
#' }
#' @export
gradient.trajectories.plot <- function(grid, f.grid, sample.points = NULL, spacing = 10) {
    if (is.null(sample.points)) {
        sample.points <- sample.points.for.visualization(nrow(f.grid), spacing)
    }

    ## Create contour plot of the function
    contour(grid$x, grid$y, f.grid, nlevels = 20)

    ## Plot gradient trajectories for sample points
    for (i in 1:nrow(sample.points)) {
        traj <- compute.gradient.trajectory(
            sample.points[i,1],
            sample.points[i,2],
            f.grid
        )

        ## Convert indices to actual coordinates
        asc.coords <- cbind(
            grid$x[traj$ascending.path[,1]],
            grid$y[traj$ascending.path[,2]]
        )
        desc.coords <- cbind(
            grid$x[traj$descending.path[,1]],
            grid$y[traj$descending.path[,2]]
        )

        ## Plot trajectories
        lines(asc.coords, col = "red")
        lines(desc.coords, col = "blue")
    }
}

morse.smale.cells1.plot <- function(grid, cell.labels) {
    ## Create color palette for different cells
    n.cells <- length(unique(as.vector(cell.labels)))
    colors <- rainbow(n.cells)

    ## Create image plot of cells
    image(grid$x, grid$y, cell.labels,
          col = colors,
          xlab = "x",
          ylab = "y",
          main = "Morse-Smale Cells")

    ## Add contour lines to show boundaries
    contour(grid$x, grid$y, cell.labels,
            add = TRUE,
            drawlabels = FALSE,
            col = "black")
}

sample.points.for.visualization <- function(axis.n.pts, spacing) {
    ## Create grid of sample points
    i.seq <- seq(1, axis.n.pts, by = spacing)
    j.seq <- seq(1, axis.n.pts, by = spacing)

    ## Create all combinations
    sample.points <- expand.grid(i = i.seq, j = j.seq)
    return(as.matrix(sample.points))
}

check.gradient.consistency <- function(f.grid) {
    issues <- list()
    is.consistent <- TRUE

    ## Check that all trajectories end at critical points
    for (i in 1:nrow(f.grid)) {
        for (j in 1:ncol(f.grid)) {
            traj <- compute.gradient.trajectory(i, j, f.grid)

            ## Check if endpoints are actually critical points
            if (!is.local.maximum(traj$local.max[1], traj$local.max[2], f.grid)) {
                issues[[length(issues) + 1]] <-
                    sprintf("Trajectory from (%d,%d) ends at non-maximum (%d,%d)",
                            i, j, traj$local.max[1], traj$local.max[2])
                is.consistent <- FALSE
            }

            if (!is.local.minimum(traj$local.min[1], traj$local.min[2], f.grid)) {
                issues[[length(issues) + 1]] <-
                    sprintf("Trajectory from (%d,%d) ends at non-minimum (%d,%d)",
                            i, j, traj$local.min[1], traj$local.min[2])
                is.consistent <- FALSE
            }
        }
    }

    return(list(is.consistent = is.consistent,
                issues = issues))
}

#' Plot Critical Points on Grid
#'
#' @description Plots local maxima and minima from grid-based critical point detection.
#'
#' @param grid Grid object from create.grid
#' @param critical.points List with 'maxima' and 'minima' matrices
#' @param max.color Color for maxima (default "red")
#' @param min.color Color for minima (default "blue")
#' @param max.pch Symbol for maxima (default 4, cross)
#' @param min.pch Symbol for minima (default 20, filled circle)
#' @param main Plot title
#' @param point.size Point size multiplier (default 2)
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param axes Logical, whether to show axes
#' @param add Logical, whether to add to existing plot
#' @param ... Additional arguments passed to plot
#' @return Invisible list with color information
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(30)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' critical <- find.critical.points(f.grid)
#' grid.critical.points.plot(grid, critical)
#' }
#' @export
grid.critical.points.plot <- function(grid, critical.points,
                                      max.color = "red",
                                      min.color = "blue",
                                      max.pch = 4,
                                      min.pch = 20,
                                      main = "Critical Points",
                                      point.size = 2,
                                      xlab = "", ylab = "",
                                      axes = FALSE,
                                      add = FALSE, ...) {
    if (!add) {
        ## Initialize plot only if not adding to existing
        plot(0, 0, type = "n",
             xlim = range(grid$x),
             ylim = range(grid$y),
             xlab = xlab, ylab = ylab,
             las = 1,
             axes = axes,
             main = main, ...)
        if (!axes) {
            box()
        }
    }

    ## Plot maxima
    if (!is.null(critical.points$maxima) && nrow(critical.points$maxima) > 0) {
        points(grid$x[critical.points$maxima[,1]],
               grid$y[critical.points$maxima[,2]],
               col = max.color, pch = max.pch, cex = point.size)
    }

    ## Plot minima
    if (!is.null(critical.points$minima) && nrow(critical.points$minima) > 0) {
        points(grid$x[critical.points$minima[,1]],
               grid$y[critical.points$minima[,2]],
               col = min.color, pch = min.pch, cex = point.size)
    }

    ## Return invisibly for potential additions to plot
    invisible(list(max.color = max.color,
                  min.color = min.color))
}

#' Plot Gradient Trajectory
#'
#' @description Plots a single gradient trajectory with ascending and descending paths.
#'
#' @param grid Grid object from create.grid
#' @param trajectory Trajectory object from compute.gradient.trajectory
#' @param ascending.color Color for ascending path (default "orange")
#' @param descending.color Color for descending path (default "purple")
#' @param line.width Line width (default 2)
#' @param add Logical, whether to add to existing plot
#' @param ... Additional arguments passed to plot
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(30)
#' f <- function(x, y) -((x-0.5)^2 + (y-0.5)^2)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' traj <- compute.gradient.trajectory(15, 15, f.grid)
#' gradient.trajectory.plot(grid, traj)
#' }
#' @export
gradient.trajectory.plot <- function(grid, trajectory,
                                   ascending.color = "orange",
                                   descending.color = "purple",
                                   line.width = 2,
                                   add = FALSE, ...) {
    ## Convert indices to coordinates
    asc.x <- grid$x[trajectory$ascending.path[,1]]
    asc.y <- grid$y[trajectory$ascending.path[,2]]
    desc.x <- grid$x[trajectory$descending.path[,1]]
    desc.y <- grid$y[trajectory$descending.path[,2]]

    if (!add) {
        ## Initialize new plot if not adding to existing
        plot(0, 0, type = "n",
             xlim = range(grid$x),
             ylim = range(grid$y),
             xlab = "x", ylab = "y",
             main = "Gradient Trajectory", ...)
    }

    ## Plot ascending path
    lines(asc.x, asc.y,
          col = ascending.color,
          lwd = line.width)

    ## Plot descending path
    lines(desc.x, desc.y,
          col = descending.color,
          lwd = line.width)

    invisible(NULL)
}

# Create a Point class using R's S3 system
create_point <- function(x, y) {
  point <- list(x = x, y = y)
  class(point) <- "Point"
  return(point)
}

#' Find Line Intersection
#'
#' @description Finds the intersection point of two line segments.
#'
#' @param P1 First point of first line segment c(x, y)
#' @param P2 Second point of first line segment c(x, y)
#' @param Q1 First point of second line segment c(x, y)
#' @param Q2 Second point of second line segment c(x, y)
#' @return Coordinates of intersection point or NULL if no intersection
#' @export
find.line.intersection <- function(P1, P2, Q1, Q2) {
    # Direction vectors
    d1 <- P2 - P1
    d2 <- Q2 - Q1

    # Cross product of direction vectors
    cross <- d1[1] * d2[2] - d1[2] * d2[1]

    # If cross product is zero, lines are parallel
    if (abs(cross) < 1e-10) {
        return(NULL)
    }

    # Solve for intersection parameters
    t <- ((Q1[1] - P1[1]) * d2[2] - (Q1[2] - P1[2]) * d2[1]) / cross

    # Calculate intersection point
    intersection <- P1 + t * d1

    return(intersection)
}

# Helper function to plot arrow on trajectory
arrow.on.trajectory.plot <- function(trajectory, index, arrow.length = 0.02,
                                   arrow.angle = 30, arrow.col = "black") {
    if(index < 2 || index >= nrow(trajectory)) return()

    # Get direction at index
    direction <- trajectory[index + 1,] - trajectory[index - 1,]
    direction <- direction / sqrt(sum(direction^2))

    # Arrow position
    arrow_pos <- trajectory[index,]

    # Calculate arrow head points
    arrow_angle_rad <- arrow.angle * pi / 180
    cos_angle <- cos(arrow_angle_rad)
    sin_angle <- sin(arrow_angle_rad)

    # Rotate direction vector for arrow heads
    left_angle <- atan2(direction[2], direction[1]) + arrow_angle_rad
    right_angle <- atan2(direction[2], direction[1]) - arrow_angle_rad

    arrow_left <- c(
        arrow_pos[1] - arrow.length * cos(left_angle),
        arrow_pos[2] - arrow.length * sin(left_angle)
    )

    arrow_right <- c(
        arrow_pos[1] - arrow.length * cos(right_angle),
        arrow_pos[2] - arrow.length * sin(right_angle)
    )

    # Draw arrow head
    lines(c(arrow_left[1], arrow_pos[1], arrow_right[1]),
          c(arrow_left[2], arrow_pos[2], arrow_right[2]),
          col = arrow.col)
}

#' Identify Destination Minimum
#'
#' @description Identifies which minimum a gradient descent trajectory converges to
#' by finding the closest known minimum to the trajectory endpoint. This function
#' is useful for classifying trajectories in Morse-Smale cell construction.
#'
#' @param trajectory Matrix of trajectory points where each row is an (x, y) coordinate
#' @param minima_coordinates Matrix of known minima locations where each row is an (x, y) coordinate
#' @param tolerance Maximum distance between trajectory endpoint and minimum to consider a match
#'
#' @return Integer index of the destination minimum (row number in minima_coordinates),
#'   or NA if no minimum is within tolerance
#'
#' @export
destination.minimum.identify <- function(trajectory, minima_coordinates, tolerance = 0.05) {
    # Get the final point of the trajectory
    final_point <- trajectory[nrow(trajectory), ]

    # Calculate distances from the final point to all minima
    distances <- apply(minima_coordinates, 1, function(min_coord) {
        sqrt(sum((final_point - min_coord)^2))
    })

    # Find the closest minimum
    min_idx <- which.min(distances)
    min_distance <- distances[min_idx]

    # Check if within tolerance
    if (min_distance <= tolerance) {
        return(min_idx)
    } else {
        warning(sprintf("Trajectory endpoint (%.3f, %.3f) is far (%.3f) from nearest minimum.",
                       final_point[1], final_point[2], min_distance))
        return(NA)
    }

}

#' Plot Function Contours
#'
#' @description Creates a contour plot of function values on a grid.
#'
#' @param grid Grid object from create.grid
#' @param f.grid Matrix of function values
#' @param main Plot title
#' @param nlevels Number of contour levels (default 20)
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param axes Logical, whether to show axes
#' @param add Logical, whether to add to existing plot
#' @param ... Additional arguments passed to contour
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(50)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' function.contours.plot(grid, f.grid)
#' }
#' @export
function.contours.plot <- function(grid, f.grid,
                                 main = "Function Contours",
                                 nlevels = 20,
                                 xlab = "", ylab = "",
                                 axes = FALSE,
                                 add = FALSE, ...) {
    if (!add) {
        ## Create new contour plot
        contour(grid$x, grid$y, f.grid,
                nlevels = nlevels,
                xlab = xlab, ylab = ylab,
                las = 1,
                axes = axes,
                main = main, ...)
    } else {
        ## Add contours to existing plot
        contour(grid$x, grid$y, f.grid,
                nlevels = nlevels,
                las = 1,
                axes = axes,
                add = TRUE, ...)
    }
    invisible(NULL)
}

#' Plot Morse Analysis
#'
#' @description Combined plot showing function contours, critical points, and a gradient trajectory.
#'
#' @param grid Grid object from create.grid
#' @param f.grid Matrix of function values
#' @param trajectory Trajectory object from compute.gradient.trajectory
#' @param critical.points List with critical points from find.critical.points
#' @param show.contours Logical, whether to show contour lines (default TRUE)
#' @param main Plot title
#' @param ... Additional arguments passed to plotting functions
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(30)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' critical <- find.critical.points(f.grid)
#' traj <- compute.gradient.trajectory(15, 20, f.grid)
#' morse.analysis.plot(grid, f.grid, traj, critical)
#' }
#' @export
morse.analysis.plot <- function(grid,
                                f.grid,
                                trajectory,
                                critical.points,
                                show.contours = TRUE,
                                main = "Morse Analysis",
                                ...) {
    ## Start with contours if requested
    if (show.contours) {
        function.contours.plot(grid, f.grid, main = main, ...)
    } else {
        ## If no contours, initialize empty plot
        plot(0, 0, type = "n",
             xlim = range(grid$x),
             ylim = range(grid$y),
             xlab = "x", ylab = "y",
             main = main, ...)
    }

    ## Add critical points
    colors <- grid.critical.points.plot(grid, critical.points, add = TRUE)

    ## Add trajectory
    gradient.trajectory.plot(grid, trajectory,
                           ascending.color = "orange",
                           descending.color = "purple",
                           add = TRUE)

    ## Add legend
    legend("topright",
           legend = c("Local Maximum", "Local Minimum",
                     "Ascending Path", "Descending Path"),
           col = c(colors$max.color, colors$min.color,
                  "orange", "purple"),
           pch = c(4, 20, NA, NA),
           lty = c(NA, NA, 1, 1),
           pt.cex = 2,
           title = "Legend")

    invisible(NULL)
}

## Function to evaluate along specific directions through a point
evaluate.along.direction <- function(f, center, direction, range = 0.4, points = 100) {
    ## Create sequence of points along the direction
    t <- seq(-range, range, length.out = points)
    x <- center[1] + direction[1] * t
    y <- center[2] + direction[2] * t

    ## Evaluate function
    values <- numeric(points)
    for(i in 1:points) {
        values[i] <- f(x[i], y[i])
    }

    return(list(t = t, values = values))
}

## ------------------------------------------------------------------
##
## using analytical gradient calculation
##
## ------------------------------------------------------------------

#' Find Ascending Trajectory from Point
#'
#' @description Traces gradient ascent path from a starting point using analytical gradient.
#'
#' @param x Starting x coordinate
#' @param y Starting y coordinate
#' @param step Step size for gradient ascent
#' @param mixture Object with gradient method
#' @param max_steps Maximum number of steps (default 1000)
#' @param grad_threshold Gradient magnitude threshold (default 1e-5)
#' @param angle_threshold Angle change threshold (default -0.5)
#' @param domain Domain boundaries as list
#' @return Matrix of (x,y) coordinates along trajectory
#' @export
find.ascending.trajectory.from.point <- function(x, y, step, mixture,
                                                max_steps = 1000, grad_threshold = 1e-5,
                                                angle_threshold = -0.5,
                                                domain = list(x_min = 0, x_max = 1,
                                                              y_min = 0, y_max = 1)) {

    # Get gradient at current point
    grad <- mixture$gradient(x, y)

    # If gradient magnitude is very small, we're at/near a critical point
    if (sqrt(sum(grad^2)) < 1e-6) return(NULL)

    # Initialize variables for finding best neighbor
    max_alignment <- -Inf
    best_idx <- NULL

    # Initialize trajectory with starting point
    trajectory <- matrix(c(x, y), nrow = 1)

    # Keep track of previous gradient for angle checking
    prev_grad <- grad / sqrt(sum(grad^2))

    for(i in 2:max_steps) {
        # Get gradient at current point
        grad <- mixture$gradient(x, y)
        grad_norm <- sqrt(sum(grad^2))

        # Stop if gradient is too small (we're at a critical point)
        if(grad_norm < grad_threshold) {
            break
        }

        # Project gradient if on boundary
        proj_grad <- project.gradient.to.boundary(x, y, grad, domain)

        # If projection killed the gradient, we're at a boundary critical point
        if(sqrt(sum(proj_grad^2)) < grad_threshold) {
            break
        }

        # Normalize gradient
        grad <- proj_grad / sqrt(sum(proj_grad^2))

        # Check angle between current and previous gradient
        angle_cos <- sum(grad * prev_grad)
        if(angle_cos < angle_threshold) {
            # Gradient direction changed significantly
            break
        }

        # Update position
        new_x <- x + step * grad[1]
        new_y <- y + step * grad[2]

        # Ensure we stay within domain
        new_x <- max(domain$x_min, min(domain$x_max, new_x))
        new_y <- max(domain$y_min, min(domain$y_max, new_y))

        # Check if we've made any progress
        if(abs(new_x - x) < 1e-10 && abs(new_y - y) < 1e-10) {
            break
        }

        x <- new_x
        y <- new_y

        # Store new position
        trajectory <- rbind(trajectory, c(x, y))

        # Update previous gradient
        prev_grad <- grad
    }

    return(trajectory)
}

#' Find Descending Trajectory from Point
#'
#' @description Traces gradient descent path from a starting point using analytical gradient.
#'
#' @param x Starting x coordinate
#' @param y Starting y coordinate
#' @param step Step size for gradient descent
#' @param mixture Object with gradient method
#' @param max_steps Maximum number of steps (default 1000)
#' @param grad_threshold Gradient magnitude threshold (default 1e-5)
#' @param angle_threshold Angle change threshold (default -0.5)
#' @param domain Domain boundaries as list
#' @return Matrix of (x,y) coordinates along trajectory
#' @export
find.descending.trajectory.from.point <- function(x, y, step = 0.01, mixture,
                                                max_steps = 1000, grad_threshold = 1e-5,
                                                angle_threshold = -0.5,
                                                domain = list(x_min = 0, x_max = 1,
                                                              y_min = 0, y_max = 1)) {
    # Create a wrapper that reverses the gradient
    descending_mixture <- list(
        f = mixture$f,
        gradient = function(x, y) -mixture$gradient(x, y)
    )

    return(find.ascending.trajectory.from.point(x, y, step, descending_mixture,
                                              max_steps, grad_threshold, angle_threshold, domain))
}

#' Create Morse-Smale Complex Object
#'
#' @description Constructs complete Morse-Smale complex from grid function values.
#'
#' @param f.grid Matrix of function values
#' @return S3 object of class "morse.smale.complex" with components:
#'   \item{cell_lookup}{List mapping cell IDs to point coordinates}
#'   \item{cell_summary}{Data frame summarizing each cell}
#'   \item{local_maxima}{Matrix of local maximum coordinates}
#'   \item{local_minima}{Matrix of local minimum coordinates}
#'   \item{find_cell}{Function to find cell ID for a grid point}
#'   \item{get_neighboring_cells}{Function to find neighboring cells}
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(30)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' complex <- create.morse.smale.complex(f.grid)
#' str(complex)
#' }
#' @export
create.morse.smale.complex <- function(f.grid) {
    n <- nrow(f.grid)

    # Find critical points
    critical.points <- find.critical.points(f.grid)
    local_maxima <- critical.points$maxima
    local_minima <- critical.points$minima

    # Create lookup table for cells
    cell_lookup <- list()

    # Process each grid point
    for (i in 1:n) {
        for (j in 1:n) {
            # Skip if this is a critical point
            if (any(apply(local_maxima, 1, function(p) all(p == c(i, j)))) ||
                any(apply(local_minima, 1, function(p) all(p == c(i, j))))) {
                next
            }

            # Compute gradient trajectory
            traj <- compute.gradient.trajectory(i, j, f.grid)

            # Create cell identifier from destination extrema
            max_idx <- which(apply(local_maxima, 1,
                                 function(p) all(p == traj$local.max)))
            min_idx <- which(apply(local_minima, 1,
                                 function(p) all(p == traj$local.min)))

            if (length(max_idx) > 0 && length(min_idx) > 0) {
                cell_id <- paste0(max_idx[1], "_", min_idx[1])

                if (is.null(cell_lookup[[cell_id]])) {
                    cell_lookup[[cell_id]] <- matrix(c(i, j), nrow = 1)
                } else {
                    cell_lookup[[cell_id]] <- rbind(cell_lookup[[cell_id]], c(i, j))
                }
            }
        }
    }

    # Create cell summary
    cell_summary <- data.frame(
        cell_id = names(cell_lookup),
        point_count = sapply(cell_lookup, nrow),
        stringsAsFactors = FALSE
    )

    # Extract max and min indices
    cell_summary$max_idx <- as.integer(sub("_.*", "", cell_summary$cell_id))
    cell_summary$min_idx <- as.integer(sub(".*_", "", cell_summary$cell_id))

    # Create Morse-Smale complex object
    morse_smale_complex <- list(
        cell_lookup = cell_lookup,
        cell_summary = cell_summary,
        local_maxima = local_maxima,
        local_minima = local_minima
    )

    # Add S3 class
    class(morse_smale_complex) <- "morse.smale.complex"

    # Add method to find which cell a point belongs to
    morse_smale_complex$find_cell <- function(i, j) {
        for (cell_id in names(cell_lookup)) {
            cell_points <- cell_lookup[[cell_id]]
            if (any(apply(cell_points, 1, function(p) all(p == c(i, j))))) {
                return(cell_id)
            }
        }
        # Check if point is a local extremum
        if (any(apply(local_maxima, 1, function(p) all(p == c(i, j))))) {
            return("local_maximum")
        }
        if (any(apply(local_minima, 1, function(p) all(p == c(i, j))))) {
            return("local_minimum")
        }
        return(NA)
    }

    # Add function to get neighboring cells
    morse_smale_complex$get_neighboring_cells <- function(cell_id) {
        # Get current cell's extrema indices
        max_idx <- as.numeric(strsplit(cell_id, "_")[[1]][1])
        min_idx <- as.numeric(strsplit(cell_id, "_")[[1]][2])

        # Find cells that share either maximum or minimum
        neighbors <- unique(c(
            # Cells sharing the same maximum
            names(cell_lookup)[grep(paste0("^", max_idx, "_"),
                                  names(cell_lookup))],
            # Cells sharing the same minimum
            names(cell_lookup)[grep(paste0("_", min_idx, "$"),
                                  names(cell_lookup))]
        ))

        neighbors <- setdiff(neighbors, cell_id)  # Remove self
        return(neighbors)
    }

    return(morse_smale_complex)
}

#' Plot Morse-Smale Cells
#'
#' @description Visualizes Morse-Smale cell decomposition with critical points.
#'
#' @param grid Grid object from create.grid
#' @param morse_smale_complex Morse-Smale complex object
#' @param f.grid Optional matrix of function values for contours
#' @param ... Additional arguments passed to plot
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(30)
#' f <- function(x, y) sin(3*x) * cos(3*y)
#' f.grid <- evaluate.function.on.grid(f, grid)
#' complex <- create.morse.smale.complex(f.grid)
#' morse.smale.cells.plot(grid, complex, f.grid)
#' }
#' @export
morse.smale.cells.plot <- function(grid, morse_smale_complex, f.grid = NULL, ...) {
    # First, create a matrix of cell labels for the entire grid
    n <- length(grid$x)
    cell.labels <- matrix(NA, n, n)

    # Fill the matrix with cell identifiers
    for (i in 1:n) {
        for (j in 1:n) {
            cell.labels[i,j] <- morse_smale_complex$find_cell(i, j)
        }
    }

    # Create a color palette that emphasizes the cell structure
    # We use distinct colors for cells and special colors for extrema
    unique_cells <- setdiff(unique(as.vector(cell.labels)),
                           c("local_maximum", "local_minimum", NA))
    n.cells <- length(unique_cells)

    # Create a harmonious color palette for cells
    cell_colors <- grDevices::hcl.colors(n.cells, "Set 2")  # Using HCL colors for better distinction

    # Create a named vector for all possible labels
    color_mapping <- c(
        setNames(cell_colors, unique_cells),
        local_maximum = "red",
        local_minimum = "blue"
    )

    # Create numeric matrix for plotting
    numeric_labels <- matrix(NA, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            label <- cell.labels[i,j]
            if (!is.na(label)) {
                numeric_labels[i,j] <- which(names(color_mapping) == label)
            }
        }
    }

    # Set up the plotting area with appropriate margins
    par(mar = c(5, 4, 4, 6))  # Adjust margins to accommodate legend

    # Create base plot
    image(grid$x, grid$y, numeric_labels,
          col = color_mapping,
          xlab = "x",
          ylab = "y",
          main = "Morse-Smale Complex", ...)

    # Add contour lines if f.grid is provided
    if (!is.null(f.grid)) {
        contour(grid$x, grid$y, f.grid,
                add = TRUE,
                drawlabels = FALSE,
                col = "black",
                lty = 2)
    }

    # Add cell boundaries
    contour(grid$x, grid$y, numeric_labels,
            add = TRUE,
            drawlabels = FALSE,
            col = "white",
            lwd = 2)

    # Mark critical points
    # Local maxima as upward-pointing triangles
    points(grid$x[morse_smale_complex$local_maxima[,1]],
           grid$y[morse_smale_complex$local_maxima[,2]],
           pch = 24,  # Up triangle
           col = "black",
           bg = "red",
           cex = 2)

    # Local minima as downward-pointing triangles
    points(grid$x[morse_smale_complex$local_minima[,1]],
           grid$y[morse_smale_complex$local_minima[,2]],
           pch = 25,  # Down triangle
           col = "black",
           bg = "blue",
           cex = 2)

    # Add informative legend
    legend_labels <- c(
        paste("Cell", seq_len(n.cells)),
        "Local Maximum",
        "Local Minimum"
    )
    legend_colors <- c(
        cell_colors,
        "red",
        "blue"
    )
    legend_pch <- c(
        rep(22, n.cells),  # Square for cells
        24,  # Up triangle for maxima
        25   # Down triangle for minima
    )

    legend("topright",
           legend = legend_labels,
           fill = c(rep(NA, n.cells), NA, NA),
           border = c(rep(NA, n.cells), NA, NA),
           pch = legend_pch,
           pt.bg = legend_colors,
           col = c(rep(NA, n.cells), "black", "black"),
           title = "Morse-Smale Components",
           inset = c(0.02, 0),
           xpd = TRUE)

    # Add cell statistics
    cell_stats <- morse_smale_complex$cell_summary
    cat("\nMorse-Smale Complex Statistics:\n")
    cat("Number of cells:", n.cells, "\n")
    cat("Points per cell:\n")
    print(data.frame(
        Cell = cell_stats$cell_id,
        Points = cell_stats$point_count,
        Max_idx = cell_stats$max_idx,
        Min_idx = cell_stats$min_idx
    ))

    invisible(NULL)
}

# Helper function to project gradient onto domain boundary
project.gradient.to.boundary <- function(x, y, grad, domain) {
    # Initialize projection flags
    at_x_min <- abs(x - domain$x_min) < 1e-10
    at_x_max <- abs(x - domain$x_max) < 1e-10
    at_y_min <- abs(y - domain$y_min) < 1e-10
    at_y_max <- abs(y - domain$y_max) < 1e-10

    # Handle corners first
    if((at_x_min && at_y_min) || (at_x_min && at_y_max) ||
       (at_x_max && at_y_min) || (at_x_max && at_y_max)) {
        return(c(0, 0))  # At corner, no movement possible
    }

    # Project gradient onto boundaries
    if(at_x_min || at_x_max) {
        grad[1] <- 0  # No x-component on vertical boundaries
    }
    if(at_y_min || at_y_max) {
        grad[2] <- 0  # No y-component on horizontal boundaries
    }

    return(grad)
}

#' Plot Critical Points
#'
#' @description Plots critical points (maxima, minima, saddles) from continuous analysis.
#'
#' @param critical_points List with maxima, minima, and saddles matrices
#' @param xlim X-axis limits
#' @param ylim Y-axis limits
#' @param main Plot title
#' @param add Logical, whether to add to existing plot
#' @param ... Additional arguments passed to plot
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' # Create a function with known critical points
#' f <- function(x, y) x^2 - y^2  # Saddle at origin
#' gradient <- function(x, y) c(2*x, -2*y)
#' grid <- create.grid(30)
#' critical <- find.critical.points.continuous(gradient, grid)
#' critical.points.plot(critical)
#' }
#' @export
critical.points.plot <- function(critical_points,
                               xlim = c(0, 1), ylim = c(0, 1),
                               main = "Critical Points Analysis",
                               add = FALSE, ...) {
    if (!add) {
        plot(NULL, xlim = xlim, ylim = ylim,
             xlab = "x", ylab = "y", main = main, ...)
    }

    # Plot maxima
    if (nrow(critical_points$maxima) > 0) {
        points(critical_points$maxima[,1], critical_points$maxima[,2],
               pch = 24, col = "black", bg = "red", cex = 2)
    }

    # Plot minima
    if (nrow(critical_points$minima) > 0) {
        points(critical_points$minima[,1], critical_points$minima[,2],
               pch = 25, col = "black", bg = "blue", cex = 2)
    }

    # Plot saddles
    if (nrow(critical_points$saddles) > 0) {
        points(critical_points$saddles[,1], critical_points$saddles[,2],
               pch = 23, col = "black", bg = "green", cex = 2)
    }

    # Add legend
    legend("topright",
           legend = c("Local Maximum", "Local Minimum", "Saddle Point"),
           pch = c(24, 25, 23),
           pt.bg = c("red", "blue", "green"),
           col = "black",
           pt.cex = 2)

    invisible(NULL)
}

#' Find Separatrices from Saddle Points
#'
#' @description Computes separatrices emanating from saddle points.
#'
#' @param saddle_point Coordinates of saddle point c(x, y)
#' @param mixture Object with gradient method
#' @param eps Small offset from saddle (default 0.01)
#' @param step Step size for trajectory (default 0.01)
#' @return List of separatrix trajectories
#' @export
find.separatrices <- function(saddle_point, mixture, eps = 0.01, step = 0.01) {
    # Function to compute trajectory from a starting point
    compute_trajectory <- function(start_point, ascending = TRUE) {
        if (ascending) {
            find.ascending.trajectory.from.point(
                start_point[1], start_point[2], step, mixture)
        } else {
            find.descending.trajectory.from.point(
                start_point[1], start_point[2], step, mixture)
        }
    }

    # Initialize separatrices list
    separatrices <- list()

    # Ascending trajectories (to maxima)
    separatrices$saddle_ascending1 <- compute_trajectory(
        saddle_point + c(eps, eps), ascending = TRUE)
    separatrices$saddle_ascending2 <- compute_trajectory(
        saddle_point + c(-eps, -eps), ascending = TRUE)

    # Descending trajectories (to minima)
    separatrices$saddle_descending1 <- compute_trajectory(
        saddle_point + c(eps, -eps), ascending = FALSE)
    separatrices$saddle_descending2 <- compute_trajectory(
        saddle_point + c(-eps, eps), ascending = FALSE)

    return(separatrices)
}

#' Plot Separatrices
#'
#' @description Plots separatrices emanating from saddle points.
#'
#' @param separatrices List of separatrix trajectories
#' @param M1_pos Position of first maximum
#' @param M2_pos Position of second maximum
#' @param line.lwd Line width (default 2)
#' @param line.lty Line type (default 2, dashed)
#' @param line.col Line color (default "black")
#' @param add Logical, whether to add to existing plot
#' @param ... Additional arguments passed to plot
#' @return Invisible NULL
#' @export
separatrices.plot <- function(separatrices, M1_pos, M2_pos,
                              line.lwd = 2,
                              line.lty = 2,
                              line.col = "black",
                              add = TRUE, ...) {
    ## Create new plot if needed
    if (!add) {
        plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             xlab = "", ylab = "", ...)
    }

    ## For boundary separatrices, connect them to the appropriate maximum
    for (edge in c("boundary_bottom", "boundary_left")) {
        if (!is.null(separatrices[[edge]])) {
            segments(x0 = separatrices[[edge]][1],
                     y0 = separatrices[[edge]][2],
                     x1 = M1_pos[1],
                     y1 = M1_pos[2],
                     lty = line.lty,
                     lwd = line.lwd,
                     col = line.col)
        }
    }

    ## Connect top and right boundary maxima to M2
    for (edge in c("boundary_top", "boundary_right")) {
        if (!is.null(separatrices[[edge]])) {
            segments(x0 = separatrices[[edge]][1],
                     y0 = separatrices[[edge]][2],
                     x1 = M2_pos[1],
                     y1 = M2_pos[2],
                     lty = line.lty,
                     lwd = line.lwd,
                     col = line.col)
        }
    }

    ## Draw line connecting M1 and M2
    segments(x0 = M1_pos[1],
             y0 = M1_pos[2],
             x1 = M2_pos[1],
             y1 = M2_pos[2],
             lty = line.lty,
             lwd = line.lwd,
             col = line.col)

    ## Draw line connecting descending trajectories
    if (!is.null(separatrices$saddle_descending1) &&
        !is.null(separatrices$saddle_descending2)) {
        d1 <- separatrices$saddle_descending1
        d1.pt <- d1[nrow(d1),]
        d2 <- separatrices$saddle_descending2
        d2.pt <- d2[nrow(d2),]
        segments(x0 = d1.pt[1],
                 y0 = d1.pt[2],
                 x1 = d2.pt[1],
                 y1 = d2.pt[2],
                 lty = line.lty,
                 lwd = line.lwd,
                 col = line.col)
    }

    invisible(NULL)
}

#' Plot Separatrices with Cells
#'
#' @description Plots Morse-Smale cells defined by separatrices using filled polygons.
#' This function visualizes the decomposition of the domain into cells based on
#' gradient flow separatrices emanating from saddle points.
#'
#' @param separatrices List of separatrix trajectories containing:
#'   \itemize{
#'     \item boundary_left: Boundary critical point on left edge
#'     \item boundary_right: Boundary critical point on right edge
#'     \item boundary_top: Boundary critical point on top edge
#'     \item boundary_bottom: Boundary critical point on bottom edge
#'     \item saddle_ascending1: First ascending trajectory from saddle
#'     \item saddle_ascending2: Second ascending trajectory from saddle
#'     \item saddle_descending1: First descending trajectory from saddle
#'     \item saddle_descending2: Second descending trajectory from saddle
#'   }
#' @param M1.pos Position of first maximum as c(x, y)
#' @param M2.pos Position of second maximum as c(x, y)
#' @param saddle_pos Position of saddle point as c(x, y)
#' @param add Logical, whether to add to existing plot (default TRUE)
#' @param cell_colors List of colors for each cell with names:
#'   \itemize{
#'     \item M1_m1: Cell connecting M1 to minimum m1
#'     \item M1_m2: Cell connecting M1 to minimum m2
#'     \item M1_m3: Cell connecting M1 to minimum m3
#'     \item M2_m2: Cell connecting M2 to minimum m2
#'     \item M2_m3: Cell connecting M2 to minimum m3
#'     \item M2_m4: Cell connecting M2 to minimum m4
#'   }
#' @param ... Additional arguments passed to plot if add = FALSE
#'
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' # Typically used with pre-computed separatrices
#' # Example with a simple two-maximum system
#' separatrices <- list(
#'   boundary_left = c(0, 0.5),
#'   boundary_right = c(1, 0.5),
#'   boundary_top = c(0.5, 1),
#'   boundary_bottom = c(0.5, 0),
#'   saddle_descending1 = matrix(c(0.5, 0.5, 0.3, 0.3), nrow = 2, byrow = TRUE),
#'   saddle_descending2 = matrix(c(0.5, 0.5, 0.7, 0.7), nrow = 2, byrow = TRUE)
#' )
#' separatrices.with.cells.plot(separatrices,
#'                             M1.pos = c(0.3, 0.7),
#'                             M2.pos = c(0.7, 0.3),
#'                             saddle_pos = c(0.5, 0.5))
#' }
#' @export
separatrices.with.cells.plot <- function(separatrices, M1.pos, M2.pos, saddle_pos,
                                       add = TRUE,
                                       cell_colors = list(
                                           "M1_m1" = rgb(0, 0.5, 0, alpha = 0.1),
                                           "M1_m2" = rgb(1, 1, 0, alpha = 0.1),
                                           "M1_m3" = rgb(1, 0, 0, alpha = 0.1),
                                           "M2_m2" = rgb(0, 0, 1, alpha = 0.1),
                                           "M2_m3" = rgb(0.6, 0, 1, alpha = 0.1),
                                           "M2_m4" = rgb(0, 1, 1, alpha = 0.1)
                                       ), ...) {

    # Create new plot if needed
    if (!add) {
        plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             xlab = "", ylab = "", ...)
    }

    ## Fill cells using polygons

    ## M1_m1 cell - connects M1 to corner minimum at (0,0)
    polygon(
        x = c(0, separatrices$boundary_left[1], M1.pos[1], separatrices$boundary_bottom[1], 0),
        y = c(0, separatrices$boundary_left[2], M1.pos[2], separatrices$boundary_bottom[2], 0),
        col = cell_colors$M1_m1,
        border = NA
    )

    ## Find saddle point as intersection of trajectories
    d1 <- separatrices$saddle_descending1
    d1.pt <- d1[nrow(d1),]
    d2 <- separatrices$saddle_descending2
    d2.pt <- d2[nrow(d2),]
    saddle.pt <- find.line.intersection(d1.pt, d2.pt, M1.pos, M2.pos)

    ## M1_m2 cell
    polygon(
        x = c(d2.pt[1], saddle.pt[1], M1.pos[1], separatrices$boundary_left[1], d2.pt[1]),
        y = c(d2.pt[2], saddle.pt[2], M1.pos[2], separatrices$boundary_left[2], d2.pt[2]),
        col = cell_colors$M1_m2,
        border = NA
    )

    ## M1_m3 cell
    polygon(
        x = c(d1.pt[1], saddle.pt[1], M1.pos[1], separatrices$boundary_bottom[1], d1.pt[1]),
        y = c(d1.pt[2], saddle.pt[2], M1.pos[2], separatrices$boundary_bottom[2], d1.pt[2]),
        col = cell_colors$M1_m3,
        border = NA
    )

    ## M2_m4 cell - connects M2 to corner minimum at (1,1)
    polygon(
        x = c(1, separatrices$boundary_top[1], M2.pos[1], separatrices$boundary_right[1], 1),
        y = c(1, separatrices$boundary_top[2], M2.pos[2], separatrices$boundary_right[2], 1),
        col = cell_colors$M2_m4,
        border = NA
    )

    ## M2_m2 cell
    polygon(
        x = c(0, separatrices$boundary_top[1], M2.pos[1], saddle.pt[1], d2.pt[1], 0),
        y = c(1, separatrices$boundary_top[2], M2.pos[2], saddle.pt[2], d2.pt[2], 1),
        col = cell_colors$M2_m2,
        border = NA
    )

    ## M2_m3 cell
    polygon(
        x = c(1, d1.pt[1], saddle.pt[1], M2.pos[1], separatrices$boundary_right[1], 1),
        y = c(0, d1.pt[2], saddle.pt[2], M2.pos[2], separatrices$boundary_right[2], 0),
        col = cell_colors$M2_m3,
        border = NA
    )

    invisible(NULL)
}

#' Plot Morse-Smale Complex from Function
#'
#' @description Creates a comprehensive visualization of the Morse-Smale complex
#' by plotting contour lines and gradient flow trajectories from a grid of sample points.
#' Ascending trajectories (to maxima) are shown in red, descending trajectories (to minima)
#' are shown in blue.
#'
#' @param grid Grid object created by create.grid containing x, y coordinates
#' @param mixture Object containing:
#'   \itemize{
#'     \item f: Function taking (x, y) and returning scalar value
#'     \item gradient: Function taking (x, y) and returning gradient vector c(dx, dy)
#'   }
#' @param n_sample_points Number of sample points along each axis for trajectory visualization
#' @param ... Additional arguments passed to function.contours.plot
#'
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' # Create a simple function with two critical points
#' grid <- create.grid(50)
#' mixture <- list(
#'   f = function(x, y) sin(3*pi*x) * cos(3*pi*y),
#'   gradient = function(x, y) {
#'     c(3*pi*cos(3*pi*x)*cos(3*pi*y),
#'       -3*pi*sin(3*pi*x)*sin(3*pi*y))
#'   }
#' )
#' morse.smale.complex.plot(grid, mixture, n_sample_points = 10)
#' }
#' @export
morse.smale.complex.plot <- function(grid, mixture, n_sample_points = 20, ...) {
    # Create base contour plot
    f_grid <- evaluate.function.on.grid(mixture$f, grid)
    function.contours.plot(grid, f_grid, main = "Morse-Smale Complex", nlevels = 10, ...)

    # Plot trajectories from sample points
    x_vals <- seq(0, 1, length.out = n_sample_points)
    y_vals <- seq(0, 1, length.out = n_sample_points)

    for(x in x_vals) {
        for(y in y_vals) {
            # Find ascending and descending trajectories
            asc_traj <- find.ascending.trajectory.from.point(x, y, 0.01, mixture)
            desc_traj <- find.descending.trajectory.from.point(x, y, 0.01, mixture)

            # Plot trajectories if they exist
            if(!is.null(asc_traj) && nrow(asc_traj) > 1) {
                lines(asc_traj[,1], asc_traj[,2], col = "red", lwd = 0.5)
            }
            if(!is.null(desc_traj) && nrow(desc_traj) > 1) {
                lines(desc_traj[,1], desc_traj[,2], col = "blue", lwd = 0.5)
            }
        }
    }

    invisible(NULL)
}

#' Plot Morse-Smale Complex from Critical Points
#'
#' @description Plots the Morse-Smale complex with emphasis on critical point structure.
#' Shows contour lines, critical points (maxima as red triangles, minima as blue triangles,
#' saddles as green diamonds), and optionally the separatrices emanating from saddle points.
#'
#' @param critical_points List containing:
#'   \itemize{
#'     \item maxima: Matrix of maximum coordinates (n x 2)
#'     \item minima: Matrix of minimum coordinates (n x 2)
#'     \item saddles: Matrix of saddle point coordinates (n x 2)
#'   }
#' @param mixture Object containing:
#'   \itemize{
#'     \item f: Function taking (x, y) and returning scalar value
#'     \item gradient: Function taking (x, y) and returning gradient vector c(dx, dy)
#'   }
#' @param grid Grid object created by create.grid for function evaluation
#' @param show_separatrices Logical, whether to compute and show separatrices from saddle points
#' @param ... Additional arguments passed to function.contours.plot
#'
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' # Example with known critical points
#' grid <- create.grid(50)
#' mixture <- list(
#'   f = function(x, y) x^2 - y^2,  # Simple saddle
#'   gradient = function(x, y) c(2*x, -2*y)
#' )
#' critical_points <- list(
#'   maxima = matrix(c(0.8, 0.2), nrow = 1),
#'   minima = matrix(c(0.2, 0.8), nrow = 1),
#'   saddles = matrix(c(0.5, 0.5), nrow = 1)
#' )
#' morse.smale.complex.plot.from.critical(critical_points, mixture, grid)
#' }
#' @export
morse.smale.complex.plot.from.critical <- function(critical_points, mixture, grid,
                                                  show_separatrices = TRUE, ...) {
    # Create base contour plot
    f_grid <- evaluate.function.on.grid(mixture$f, grid)
    function.contours.plot(grid, f_grid, main = "Morse-Smale Complex", ...)

    # Add critical points
    critical.points.plot(critical_points, add = TRUE)

    # Add separatrices from saddle points
    if(show_separatrices && nrow(critical_points$saddles) > 0) {
        for(i in 1:nrow(critical_points$saddles)) {
            saddle <- critical_points$saddles[i,]
            seps <- find.separatrices(saddle, mixture)

            # Plot each separatrix
            for(sep_name in names(seps)) {
                if(!is.null(seps[[sep_name]]) && nrow(seps[[sep_name]]) > 1) {
                    sep <- seps[[sep_name]]
                    if(grepl("ascending", sep_name)) {
                        lines(sep[,1], sep[,2], col = "darkred", lwd = 2)
                    } else {
                        lines(sep[,1], sep[,2], col = "darkblue", lwd = 2)
                    }
                }
            }
        }
    }

    invisible(NULL)
}


#' Plot Morse-Smale Trajectories from Point
#'
#' @description Visualizes gradient flow trajectories (both ascending and descending)
#' from a specific starting point. Shows the trajectories on a contour plot background,
#' marks the starting point and endpoints (critical points), and optionally adds
#' directional arrows along the trajectories.
#'
#' @param x Starting x coordinate
#' @param y Starting y coordinate
#' @param mixture Object containing:
#'   \itemize{
#'     \item f: Function taking (x, y) and returning scalar value
#'     \item gradient: Function taking (x, y) and returning gradient vector c(dx, dy)
#'   }
#' @param grid Grid object created by create.grid for contour plot background
#' @param step Step size for gradient flow integration (default 0.01)
#' @param show.arrows Logical, whether to show directional arrows on trajectories
#' @param use.custom.legend Logical, whether to use custom legend.
#' @param ... Additional arguments passed to function.contours.plot
#'
#' @return Invisible list containing:
#'   \item{ascending}{Matrix of points along ascending trajectory}
#'   \item{descending}{Matrix of points along descending trajectory}
#'
#' @examples
#' \dontrun{
#' grid <- create.grid(50)
#' mixture <- list(
#'   f = function(x, y) -((x-0.3)^2 + (y-0.7)^2) - 2*((x-0.7)^2 + (y-0.3)^2),
#'   gradient = function(x, y) {
#'     c(-2*(x-0.3) - 4*(x-0.7), -2*(y-0.7) - 4*(y-0.3))
#'   }
#' )
#' # Plot trajectories from a point
#' result <- plot.morse.smale.trajectories.from.point(0.5, 0.5, mixture, grid)
#' }
#' @export
plot.morse.smale.trajectories.from.point <- function(x,
                                                     y,
                                                     mixture,
                                                     grid,
                                                     step = 0.01,
                                                     show.arrows = TRUE,
                                                     use.custom.legend = FALSE,
                                                    ...) {
    # Create base contour plot
    f_grid <- evaluate.function.on.grid(mixture$f, grid)
    function.contours.plot(grid, f_grid,
                         main = sprintf("Trajectories from (%.2f, %.2f)", x, y), ...)

    # Mark starting point
    points(x, y, pch = 19, col = "black", cex = 1.5)

    # Compute and plot ascending trajectory
    asc_traj <- find.ascending.trajectory.from.point(x, y, step, mixture)
    if(!is.null(asc_traj) && nrow(asc_traj) > 1) {
        lines(asc_traj[,1], asc_traj[,2], col = "red", lwd = 2)
        if(show.arrows && nrow(asc_traj) > 5) {
            mid_idx <- floor(nrow(asc_traj)/2)
            arrow.on.trajectory.plot(asc_traj, mid_idx, arrow.col = "red")
        }
        # Mark endpoint
        end_pt <- asc_traj[nrow(asc_traj),]
        points(end_pt[1], end_pt[2], pch = 24, col = "black", bg = "red", cex = 1.5)
    }

    # Compute and plot descending trajectory
    desc_traj <- find.descending.trajectory.from.point(x, y, step, mixture)
    if(!is.null(desc_traj) && nrow(desc_traj) > 1) {
        lines(desc_traj[,1], desc_traj[,2], col = "blue", lwd = 2)
        if(show.arrows && nrow(desc_traj) > 5) {
            mid_idx <- floor(nrow(desc_traj)/2)
            arrow.on.trajectory.plot(desc_traj, mid_idx, arrow.col = "blue")
        }
        # Mark endpoint
        end_pt <- desc_traj[nrow(desc_traj),]
        points(end_pt[1], end_pt[2], pch = 25, col = "black", bg = "blue", cex = 1.5)
    }

    if (use.custom.legend) {
        ## custom legend
        leg_x <- par("usr")[2] * 0.75  # Adjust as needed
        leg_y <- par("usr")[4] * 0.95  # Adjust as needed
        leg_spacing <- diff(par("usr")[3:4]) * 0.05

                                        # Background box (optional)
        rect(leg_x - 0.05, leg_y - leg_spacing * 3.5,
             leg_x + 0.25, leg_y + leg_spacing * 0.5,
             col = "white", border = "black")

                                        # Start point
        points(leg_x, leg_y, pch = 19, col = "black", cex = 1.2)
        text(leg_x + 0.03, leg_y, "Start", adj = 0, cex = 0.9)

                                        # Ascent
        leg_y <- leg_y - leg_spacing
        lines(c(leg_x - 0.02, leg_x + 0.02), c(leg_y, leg_y), col = "red", lwd = 2)
        points(leg_x, leg_y, pch = 24, col = "black", bg = "red", cex = 1.2)
        text(leg_x + 0.03, leg_y, "Ascent to Max", adj = 0, cex = 0.9)

                                        # Descent
        leg_y <- leg_y - leg_spacing
        lines(c(leg_x - 0.02, leg_x + 0.02), c(leg_y, leg_y), col = "blue", lwd = 2)
        points(leg_x, leg_y, pch = 25, col = "black", bg = "blue", cex = 1.2)
        text(leg_x + 0.03, leg_y, "Descent to Min", adj = 0, cex = 0.9)
    } else {
        ## Add legend with colored lines
        legend("topright",
               legend = c("Start", "Ascent to Max", "Descent to Min"),
               pch = c(19, NA, NA),                 # Only show point for "Start"
               col = c("black", "red", "blue"),     # Colors for points and lines
               pt.bg = c("black", NA, NA),          # Fill only for "Start"
               lty = c(NA, 1, 1),                   # Lines for trajectories
               lwd = c(NA, 2, 2),                   # Line widths
               merge = FALSE)                       # Don't merge point and line

        ## Add the endpoint symbols separately
        legend("topright",
               legend = c("", "", ""),              # Empty labels
               pch = c(NA, 24, 25),                 # Triangle symbols
               col = c(NA, "black", "black"),       # Black borders
               pt.bg = c(NA, "red", "blue"),        # Colored fills
               lty = 0,                             # No lines
               bty = "n",                           # No box
               inset = c(0.12, 0))                  # Adjust position to align
    }

    invisible(list(ascending = asc_traj, descending = desc_traj))
}

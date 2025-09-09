#' Morse-Smale Complex Analysis for 2D Smooth Functions
#'
#' @description
#' Functions for computing and visualizing Morse-Smale complexes for smooth scalar functions
#' defined over \eqn{\[0,1\]^2} or other 2D rectangular domains. The implementation uses a discrete
#' gradient flow approximation on a uniform grid to construct Morse-Smale cells.
#'
#' The main functionalities include:
#' * Computing gradient trajectories
#' * Identifying critical points (local minima and maxima)
#' * Constructing Morse-Smale cells
#' * Visualizing gradient flows and cell decompositions
#'
#' @details
#' The gradient flow is approximated on a discrete uniform grid. For each grid point,
#' the algorithm considers its immediate neighbors (up to 4 for interior points,
#' 3 for edge points, 2 for corner points) to determine the direction of steepest
#' ascent/descent.
#'
#' Morse-Smale cells are constructed by grouping points whose gradient trajectories
#' flow to the same (minimum, maximum) pair of critical points.
#'
#' @section Warning:
#' The accuracy of the Morse-Smale complex approximation depends on the grid resolution
#' and the smoothness of the input function. Higher grid resolutions provide better
#' approximations but increase computation time.
#'
#' @author Pawel Gajer
#'
#' @references
#' Milnor, J. (1963). Morse Theory. Princeton University Press.
#'
#' Edelsbrunner, H., & Harer, J. (2010). Computational Topology: An Introduction.
#' American Mathematical Society.
#' # A more accessible introduction with biological applications
#'
#' Singh, G., Mémoli, F., & Carlsson, G. (2007). Topological Methods for the
#' Analysis of High Dimensional Data Sets and 3D Object Recognition. Eurographics
#' Symposium on Point-Based Graphics, 91-100.
#' # Introduces Mapper algorithm, widely used in biological data analysis
#'
#' Nicolau, M., Levine, A. J., & Carlsson, G. (2011). Topology based data
#' analysis identifies a subgroup of breast cancers with a unique mutational
#' profile and excellent survival. Proceedings of the National Academy of
#' Sciences, 108(17), 7265-7270.
#' # Classic example of applying topological methods to cancer genomics
#'
#' Yao, Y., Sun, J., Huang, X., Bowman, G. R., Singh, G., Lesnick, M., Pande, V. S.,
#' Carlsson, G. (2009). Topological methods for exploring low-density states in
#' biomolecular folding pathways. The Journal of Chemical Physics, 130(14), 144115.
#' # Application to protein folding landscapes
#'
#' Cang, Z., & Wei, G. W. (2017). TopologyNet: Topology based deep convolutional
#' neural networks for biomolecular property predictions. PLOS Computational Biology,
#' 13(7), e1005690.
#' # Modern application combining topology with machine learning
#'
#' Patania, A., Vaccarino, F., & Petri, G. (2017). Topological analysis of data.
#' EPJ Data Science, 6(1), 7.
#' # Accessible review of topological data analysis methods
#'
#' These references progress from theoretical foundations (Milnor) through computational approaches (Edelsbrunner & Harer) to specific biological applications. Each provides different insights:
#'
#' 1. Edelsbrunner & Harer provides a more computational perspective and serves as a bridge between pure theory and applications
#' 2. The Singh et al. paper introduces the Mapper algorithm, which has become widely used in biological data analysis
#' 3. The Nicolau et al. paper demonstrates a successful application in cancer research
#' 4. The Yao et al. paper shows application to protein dynamics
#' 5. The Cang & Wei paper shows how these methods can be integrated with modern machine learning
#' 6. The Patania et al. paper provides an accessible overview of the field
#'
#' @examples
#' # Create a test function and grid
#' axis.n.pts <- 50
#' grid <- create_grid(axis.n.pts)
#' f <- function(x, y) -4 * ((x - 0.5)^2 + (y - 0.5)^2) + 1
#' f.grid <- evaluate_function_on_grid(f, grid)
#'
#' # Compute and visualize Morse-Smale cells
#' cell_labels <- compute_morse_smale_cells(f.grid)
#' plot_morse_smale_cells(grid, cell_labels)
#'
#' @importFrom graphics contour image lines par title
#' @importFrom grDevices rainbow

## Core data structure functions
create.grid <- function(axis.n.pts) {
    ## Creates uniform 2D grid over [0,1]^2
    ## Returns: list containing x, y vectors and grid data.frame
    x <- seq(0, 1, length.out = axis.n.pts)
    y <- seq(0, 1, length.out = axis.n.pts)
    grid <- expand.grid(x = x, y = y)
    return(list(x = x, y = y, grid = grid))
}


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

evaluate.function.on.grid <- function(f, grid) {
    # Initialize output matrix
    f.grid <- matrix(0, nrow = length(grid$x), ncol = length(grid$y))

    # Evaluate function at each grid point
    for(i in 1:length(grid$x)) {
        for(j in 1:length(grid$y)) {
            f.grid[i,j] <- f(grid$x[i], grid$y[j])
        }
    }

    return(f.grid)
}

evaluate.mixture.on.grid <- function(mixture, grid) {
    ## Evaluates function f on the grid
    ## Returns: matrix of function values
    f.grid <- matrix(mixture$f(grid$grid$x, grid$grid$y),
                    nrow = length(grid$x))
    return(f.grid)
}

## Neighbor analysis functions
get.neighbors <- function(i, j, axis.n.pts) {
    ## Initialize empty list for neighbors
    neighbors <- list()

    ## Check all possible neighbors while respecting boundaries
    if (i > 1) neighbors[[length(neighbors) + 1]] <- c(i-1, j)      ## Top
    if (i < axis.n.pts) neighbors[[length(neighbors) + 1]] <- c(i+1, j)  ## Bottom
    if (j > 1) neighbors[[length(neighbors) + 1]] <- c(i, j-1)      ## Left
    if (j < axis.n.pts) neighbors[[length(neighbors) + 1]] <- c(i, j+1)  ## Right

    ## Convert list to matrix
    return(do.call(rbind, neighbors))
}

# Helper function to check critical points
check.critical.point.values <- function(f.grid) {
    # Print values around the center point where we expect a maximum
    mid_i <- nrow(f.grid) %/% 2
    mid_j <- ncol(f.grid) %/% 2

    cat("Center point value:", f.grid[mid_i, mid_j], "\n")
    cat("Neighboring values:\n")
    cat("Top:", f.grid[mid_i-1, mid_j], "\n")
    cat("Bottom:", f.grid[mid_i+1, mid_j], "\n")
    cat("Left:", f.grid[mid_i, mid_j-1], "\n")
    cat("Right:", f.grid[mid_i, mid_j+1], "\n")
}

check.point.values <- function(i, j, f.grid) {
    cat("Point: (", i, ",", j, ")\n")
    cat("Value: ", f.grid[i, j], "\n")
    cat("Neighboring values:\n")
    cat("Top:", f.grid[i-1, j], "\n")
    cat("Bottom:", f.grid[i+1, j], "\n")
    cat("Left:", f.grid[i, j-1], "\n")
    cat("Right:", f.grid[i, j+1], "\n")
}


#' Find Grid Indices of Critical Points
#'
#' This function identifies the indices in a grid where critical points occur by checking
#' each point against its neighbors for local maxima and minima conditions.
#'
#' @param f.grid Matrix of function values evaluated on a grid
#' @return List containing matrices of indices for maxima and minima
find.indices.of.grid.critical.points <- function(f.grid) {
    maxima <- matrix(ncol = 2, nrow = 0)
    minima <- matrix(ncol = 2, nrow = 0)

    for (i in 1:nrow(f.grid)) {
        for (j in 1:ncol(f.grid)) {
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

#' Find Critical Points in Grid Coordinates
#'
#' This function converts grid indices to actual coordinates and performs additional
#' analysis of the critical points.
#'
#' @param grid List containing x and y vectors defining the coordinate grid
#' @param f.grid Matrix of function values evaluated on the grid
#' @return List containing coordinates of critical points and their properties
find.grid.critical.points <- function(grid, f.grid) {
    # First find indices
    indices <- find.indices.of.grid.critical.points(f.grid)

    # Convert indices to coordinates
    maxima <- matrix(ncol = 2, nrow = nrow(indices$maxima))
    minima <- matrix(ncol = 2, nrow = nrow(indices$minima))

    # Convert maxima indices to coordinates
    for(i in 1:nrow(indices$maxima)) {
        maxima[i,] <- c(
            grid$x[indices$maxima[i,1]],
            grid$y[indices$maxima[i,2]]
        )
    }

    # Convert minima indices to coordinates
    for(i in 1:nrow(indices$minima)) {
        minima[i,] <- c(
            grid$x[indices$minima[i,1]],
            grid$y[indices$minima[i,2]]
        )
    }

    # Return both coordinate and index information
    return(list(
        coordinates = list(maxima = maxima, minima = minima),
        indices = indices
    ))
}

# Modified local maximum detection
is.local.maximum <- function(i, j, f.grid) {
    neighbors <- get.neighbors(i, j, nrow(f.grid))
    center.value <- f.grid[i, j]

    # Check if center value is greater than or equal to all neighbors
    # and strictly greater than at least one neighbor
    all.leq <- TRUE
    any.strictly.less <- FALSE

    for (k in 1:nrow(neighbors)) {
        neighbor.value <- f.grid[neighbors[k,1], neighbors[k,2]]
        if (neighbor.value > center.value) {
            all.leq <- FALSE
            break
        }
        if (neighbor.value < center.value) {
            any.strictly.less <- TRUE
        }
    }

    return(all.leq && any.strictly.less)
}

# Modified local minimum detection
is.local.minimum <- function(i, j, f.grid) {
    neighbors <- get.neighbors(i, j, nrow(f.grid))
    center.value <- f.grid[i, j]

    # Check if center value is less than or equal to all neighbors
    # and strictly less than at least one neighbor
    all.geq <- TRUE
    any.strictly.greater <- FALSE

    for (k in 1:nrow(neighbors)) {
        neighbor.value <- f.grid[neighbors[k,1], neighbors[k,2]]
        if (neighbor.value < center.value) {
            all.geq <- FALSE
            break
        }
        if (neighbor.value > center.value) {
            any.strictly.greater <- TRUE
        }
    }

    return(all.geq && any.strictly.greater)
}

## Gradient trajectory functions
find.next.ascending.point1 <- function(i, j, f.grid) {
    neighbors <- get.neighbors(i, j, nrow(f.grid))
    center.value <- f.grid[i, j]

    ## Find neighbor with maximum value
    max.value <- center.value
    max.idx <- NULL

    for (k in 1:nrow(neighbors)) {
        neighbor.value <- f.grid[neighbors[k,1], neighbors[k,2]]
        if (neighbor.value > max.value) {
            max.value <- neighbor.value
            max.idx <- k
        }
    }

    ## Return NULL if no ascending direction (local maximum)
    if (is.null(max.idx)) return(NULL)
    return(neighbors[max.idx,])
}

find.next.descending.point1 <- function(i, j, f.grid) {
    neighbors <- get.neighbors(i, j, nrow(f.grid))
    center.value <- f.grid[i, j]

    ## Find neighbor with minimum value
    min.value <- center.value
    min.idx <- NULL

    for (k in 1:nrow(neighbors)) {
        neighbor.value <- f.grid[neighbors[k,1], neighbors[k,2]]
        if (neighbor.value < min.value) {
            min.value <- neighbor.value
            min.idx <- k
        }
    }

    ## Return NULL if no descending direction (local minimum)
    if (is.null(min.idx)) return(NULL)
    return(neighbors[min.idx,])
}

compute.gradient.trajectory1 <- function(i, j, f.grid) {
    ## Initialize paths
    ascending.path <- matrix(c(i, j), ncol = 2)
    descending.path <- matrix(c(i, j), ncol = 2)

    ## Compute ascending path
    current <- c(i, j)
    while (!is.null(current)) {
        next.point <- find.next.ascending.point(current[1], current[2], f.grid)
        if (!is.null(next.point)) {
            ascending.path <- rbind(ascending.path, next.point)
            current <- next.point
        } else {
            break
        }
    }
    local.max <- current

    ## Compute descending path
    current <- c(i, j)
    while (!is.null(current)) {
        next.point <- find.next.descending.point(current[1], current[2], f.grid)
        if (!is.null(next.point)) {
            descending.path <- rbind(descending.path, next.point)
            current <- next.point
        } else {
            break
        }
    }
    local.min <- current

    return(list(
        ascending.path = ascending.path,
        descending.path = descending.path,
        local.max = local.max,
        local.min = local.min
    ))
}

#' Find and Classify Critical Points of a 2D Gaussian Mixture
#'
#' @description
#' This function identifies and classifies all critical points (maxima, minima, and saddle points)
#' of a two-dimensional Gaussian mixture function. Critical points are locations where the gradient
#' vanishes, and their classification is determined by analyzing the eigenvalues of the Hessian matrix.
#'
#' @details
#' The function uses three steps to identify and classify critical points:
#'
#' 1. Critical Point Detection:
#'    - Scans the grid for points where the gradient norm is below the threshold
#'    - Uses numerical approximations to compute gradient values
#'
#' 2. Hessian Analysis:
#'    - Computes the Hessian matrix using numerical second derivatives
#'    - Calculates eigenvalues using the characteristic equation
#'
#' 3. Classification:
#'    - Maxima: Both eigenvalues negative (surface curves downward in all directions)
#'    - Minima: Both eigenvalues positive (surface curves upward in all directions)
#'    - Saddle Points: Eigenvalues have opposite signs (surface curves up in one direction, down in another)
#'
#' The eigenvalues themselves provide additional information:
#' - Their magnitude indicates the "steepness" of the surface in principal directions
#' - Their ratio indicates how "symmetric" the critical point is
#' - For saddle points, they show which direction is ascending vs descending
#'
#' @param grid List containing x and y vectors defining the search grid
#' @param mixture List containing the function and its gradient
#' @param threshold Numeric value for gradient norm threshold (default: 1e-5)
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item maxima: Matrix of (x,y) coordinates for local maxima
#'   \item minima: Matrix of (x,y) coordinates for local minima
#'   \item saddles: Matrix of (x,y) coordinates for saddle points
#'     with attribute "eigenvalues" containing corresponding eigenvalue pairs
#' }
#'
#' @examples
#' # Create a mixture with three Gaussians that will generate a saddle point
#' mixture <- create.gaussian.mixture(
#'   x1 = 0.3, y1 = 0.3, A1 = 1.0,
#'   x2 = 0.7, y2 = 0.7, A2 = 1.0,
#'   x3 = 0.5, y3 = 0.5, A3 = -0.5  # Negative amplitude creates saddle
#' )
#'
#' grid <- create.grid(100)
#' critical_points <- find.critical.points(grid, mixture)
#'
#' # Analyze saddle point characteristics
#' eigenvals <- attr(critical_points$saddles, "eigenvalues")
#' print(paste("Steepest ascent direction magnitude:", max(abs(eigenvals))))
debugging.find.critical.points.with.gradient <- function(grid, mixture, threshold = 1e-5) {
    # Initialize empty matrices for storing all types of critical points
    maxima <- matrix(ncol = 2, nrow = 0)
    minima <- matrix(ncol = 2, nrow = 0)
    saddles <- matrix(ncol = 2, nrow = 0)

    # Add debug counters
    points_checked <- 0
    small_gradients <- 0
    grad.norms <- c()

    # For each grid point
    for (i in 2:(length(grid$x) - 1)) {
        for (j in 2:(length(grid$y) - 1)) {
            points_checked <- points_checked + 1

            # Convert grid indices to actual coordinates
            x <- grid$x[i]
            y <- grid$y[j]

            # Calculate gradient at this point
            grad <- mixture$gradient(x, y)
            grad_norm <- sqrt(sum(grad^2))
            grad.norms <- c(grad.norms, grad_norm)

            # Debug print every 1000 points
            if (points_checked %% 1000 == 0) {
                cat(sprintf("Checking point %d: (%.3f, %.3f) - gradient norm: %.6f\n",
                          points_checked, x, y, grad_norm))
            }

            # If gradient magnitude is small enough, we have a critical point
            if (grad_norm < threshold) {
                small_gradients <- small_gradients + 1
                cat(sprintf("Found small gradient at (%.3f, %.3f) - norm: %.6f\n",
                          x, y, grad_norm))

                # Compute Hessian matrix elements
                h11 <- numerical.second.derivative(mixture$gradient, x, y, 1, 1, 1e-5)
                h12 <- numerical.second.derivative(mixture$gradient, x, y, 1, 2, 1e-5)
                h22 <- numerical.second.derivative(mixture$gradient, x, y, 2, 2, 1e-5)

                cat(sprintf("Hessian at (%.3f, %.3f):\n[%.6f %.6f]\n[%.6f %.6f]\n",
                          x, y, h11, h12, h12, h22))

                # Compute quantities needed for eigenvalue calculation
                trace <- h11 + h22
                det <- h11 * h22 - h12 * h12
                discriminant <- trace * trace - 4 * det

                if (discriminant >= 0) {
                    # Calculate eigenvalues
                    lambda1 <- (trace + sqrt(discriminant)) / 2
                    lambda2 <- (trace - sqrt(discriminant)) / 2

                    cat(sprintf("Eigenvalues: %.6f, %.6f\n", lambda1, lambda2))

                    # Classify based on eigenvalue signs
                    if (lambda1 < 0 && lambda2 < 0) {
                        maxima <- rbind(maxima, c(x, y))
                        cat("Classified as maximum\n")
                    } else if (lambda1 > 0 && lambda2 > 0) {
                        minima <- rbind(minima, c(x, y))
                        cat("Classified as minimum\n")
                    } else if (lambda1 * lambda2 < 0) {
                        saddles <- rbind(saddles, c(x, y))
                        cat("Classified as saddle point\n")
                    }
                }
            }
        }
    }

    # Print summary statistics
    cat(sprintf("\nSummary:\n"))
    cat(sprintf("Total points checked: %d\n", points_checked))
    cat(sprintf("Points with small gradients: %d\n", small_gradients))
    cat(sprintf("Maxima found: %d\n", nrow(maxima)))
    cat(sprintf("Minima found: %d\n", nrow(minima)))
    cat(sprintf("Saddle points found: %d\n", nrow(saddles)))

    return(list(
        maxima = maxima,
        minima = minima,
        saddles = saddles,
        grad.norms = grad.norms
    ))
}



find.critical.points.with.gradient <- function(grid, mixture, threshold = 1e-5) {
    # Initialize empty matrices for storing all types of critical points
    maxima <- matrix(ncol = 2, nrow = 0)
    minima <- matrix(ncol = 2, nrow = 0)
    saddles <- matrix(ncol = 2, nrow = 0)  # New matrix for saddle points

    # For each grid point
    for (i in 2:(length(grid$x) - 1)) {
        for (j in 2:(length(grid$y) - 1)) {
            # Convert grid indices to actual coordinates
            x <- grid$x[i]
            y <- grid$y[j]

            # Calculate gradient at this point
            grad <- mixture$gradient(x, y)
            grad_norm <- sqrt(sum(grad^2))

            # If gradient magnitude is small enough, we have a critical point
            if (grad_norm < threshold) {
                # Compute Hessian matrix elements
                h11 <- numerical.second.derivative(mixture$gradient, x, y, 1, 1, 1e-5)
                h12 <- numerical.second.derivative(mixture$gradient, x, y, 1, 2, 1e-5)
                h22 <- numerical.second.derivative(mixture$gradient, x, y, 2, 2, 1e-5)

                # Compute quantities needed for eigenvalue calculation
                trace <- h11 + h22
                det <- h11 * h22 - h12 * h12
                discriminant <- trace * trace - 4 * det

                if (discriminant >= 0) {
                    # Calculate eigenvalues using quadratic formula
                    lambda1 <- (trace + sqrt(discriminant)) / 2
                    lambda2 <- (trace - sqrt(discriminant)) / 2

                    # Classify based on eigenvalue signs
                    if (lambda1 < 0 && lambda2 < 0) {
                        maxima <- rbind(maxima, c(x, y))
                    } else if (lambda1 > 0 && lambda2 > 0) {
                        minima <- rbind(minima, c(x, y))
                    } else if (lambda1 * lambda2 < 0) {
                        # If eigenvalues have opposite signs (their product is negative)
                        saddles <- rbind(saddles, c(x, y))
                    }

                    # Store eigenvalues with each point for additional analysis
                    attr(saddles, "eigenvalues") <- cbind(
                        attr(saddles, "eigenvalues"),
                        c(lambda1, lambda2)
                    )
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
        # ∂²f/∂x∂y
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
plot.gradient.trajectories <- function(grid, f.grid, sample.points = NULL, spacing = 10) {
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

plot.morse.smale.cells1 <- function(grid, cell.labels) {
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

    return(list(
        is.consistent = is.consistent,
        issues = issues
    ))
}

##
## Plotting utilities for critical points and trajectories
##

# Plot with different colors for each type
plot.critical.points <- function(critical.points,
                                 add = FALSE,
                                 xlab = "", ylab = "",
                                 max.pch = 24,
                                 min.pch = 25,
                                 saddle.pch = 18,
                                 max.col = "red",
                                 min.col = "blue",
                                 saddle.col = "purple",
                                 point.size = 2,
                                 with.box = TRUE,
                                 axes = FALSE) {
    if (!add) {
        # Set up the plot
        plot(0, 0, type = "n", xlim = c(0,1), ylim = c(0,1), las = 1, axes = axes, xlab = xlab, ylab = ylab)
        if (!axes && with.box) {
            box()
        }
    }

    # Plot each type of critical point with different symbols/colors
    if (nrow(critical.points$coordinates$maxima) > 0) {
        points(critical.points$coordinates$maxima,
               col = max.col,
               bg  = max.col,
               cex = point.size,
               pch = max.pch)  # Filled circles for maxima
    }
    if (nrow(critical.points$coordinates$minima) > 0) {
        points(critical.points$coordinates$minima,
               col = min.col,
               bg  = min.col,
               cex = point.size,
               pch = min.pch)  # Filled circles for minima
    }
    if ("saddles" %in% names(critical.points$coordinates) && nrow(critical.points$coordinates$saddles) > 0) {
        points(critical.points$coordinates$saddles,
               col = saddle.col,
               bg  = saddle.col,
               cex = point.size,
               pch = saddle.pch)  # Diamond for saddles
    }
}

#' Draw an arrow at the midpoint of a line segment
#'
#' This function draws an arrow at the exact middle of a line segment connecting
#' two points (x0,y0) and (x1,y1). The arrow is drawn in the direction from the
#' start point to the end point.
#'
#' @param x0 Numeric value for the x-coordinate of the starting point
#' @param y0 Numeric value for the y-coordinate of the starting point
#' @param x1 Numeric value for the x-coordinate of the ending point
#' @param y1 Numeric value for the y-coordinate of the ending point
#' @param eps Numeric value controlling the length of the arrow relative to the total line length (default: 0.05)
#' @param col Color for the arrow (default: "gray40")
#' @param length Size of the arrow head (default: 0.1)
#'
#' @return No return value, called for side effect of drawing an arrow
#'
#' @examples
#' # Draw a mid-arrow from (0,0) to (1,1)
#' plot(c(0,1), c(0,1), type = "n", xlab = "x", ylab = "y")
#' draw.mid.arrow(0, 0, 1, 1)
#'
#' @export
draw.mid.arrow <- function(x0, y0, x1, y1, eps = 0.05, col = "gray40", length = 0.1) {
    midpoint_x <- (x0 +x1)/2
    midpoint_y <- (y0 + y1)/2
    ## Calculate direction vectors along each segment
    dir_x <-x1 - x0
    dir_y <- y1 - y0
    ## Calculate length of each segment for normalization
    segment_lengths <- sqrt(dir_x^2 + dir_y^2)
    ## Create unit direction vectors
    unit_dir_x <- dir_x / segment_lengths
    unit_dir_y <- dir_y / segment_lengths
    ## Calculate points eps before and after midpoints along the line direction
    before_x <- midpoint_x - eps * unit_dir_x * segment_lengths
    before_y <- midpoint_y - eps * unit_dir_y * segment_lengths
    after_x <- midpoint_x + eps * unit_dir_x * segment_lengths
    after_y <- midpoint_y + eps * unit_dir_y * segment_lengths
    ## Draw arrows from points before midpoints to points after midpoints
    arrows(before_x, before_y, after_x, after_y, length = length, col = col)
}

#' Draw an arrow at any specified fraction of a line segment
#'
#' This function draws an arrow at a point that is a fraction 'p' of the distance
#' from the starting point (x0,y0) to the ending point (x1,y1). The arrow is drawn
#' in the direction from the start point to the end point.
#'
#' @param x0 Numeric value for the x-coordinate of the starting point
#' @param y0 Numeric value for the y-coordinate of the starting point
#' @param x1 Numeric value for the x-coordinate of the ending point
#' @param y1 Numeric value for the y-coordinate of the ending point
#' @param p Numeric value between 0 and 1 indicating the fraction of the distance
#'   from (x0,y0) to (x1,y1) where the arrow should be placed (default: 0.5, which
#'   places it at the midpoint)
#' @param eps Numeric value controlling the length of the arrow relative to the total line length (default: 0.05)
#' @param col Color for the arrow (default: "gray40")
#' @param length Size of the arrow head (default: 0.1)
#'
#' @return No return value, called for side effect of drawing an arrow
#'
#' @examples
#' # Draw arrows at different positions along a line
#' plot(c(0,1), c(0,1), type = "n", xlab = "x", ylab = "y")
#' segments(0, 0, 1, 1, lty = 2) # Draw the full line segment as reference
#' draw.p.mid.arrow(0, 0, 1, 1, p = 0.25, col = "red")  # Arrow at 25% of the way
#' draw.p.mid.arrow(0, 0, 1, 1, p = 0.5, col = "blue")  # Arrow at midpoint
#' draw.p.mid.arrow(0, 0, 1, 1, p = 0.75, col = "green") # Arrow at 75% of the way
#'
#' @export
draw.p.mid.arrow <- function(x0, y0, x1, y1, p = 0.5, eps = 0.05, col = "gray40", length = 0.1) {
    ## Calculate the point that is p fraction of the way from (x0,y0) to (x1,y1)
    midpoint_x <- x0 + p * (x1 - x0)
    midpoint_y <- y0 + p * (y1 - y0)

    ## Calculate direction vectors along each segment
    dir_x <- x1 - x0
    dir_y <- y1 - y0

    ## Calculate length of each segment for normalization
    segment_lengths <- sqrt(dir_x^2 + dir_y^2)

    ## Create unit direction vectors
    unit_dir_x <- dir_x / segment_lengths
    unit_dir_y <- dir_y / segment_lengths

    ## Calculate points eps before and after midpoints along the line direction
    before_x <- midpoint_x - eps * unit_dir_x * segment_lengths
    before_y <- midpoint_y - eps * unit_dir_y * segment_lengths
    after_x <- midpoint_x + eps * unit_dir_x * segment_lengths
    after_y <- midpoint_y + eps * unit_dir_y * segment_lengths

    ## Draw arrows from points before midpoints to points after midpoints
    arrows(before_x, before_y, after_x, after_y, length = length, col = col)
}

plot.grid.critical.points <- function(grid, critical.points,
                               main = "",
                               max.color = "red",
                               min.color = "blue",
                               max.pch = 19,   # x
                               min.pch = 19,   # solid dot
                               point.size = 2,
                               xlab = "", ylab = "",
                               axes = FALSE,
                               add = FALSE) {  # Added add parameter
    if (!add) {
        ## Initialize plot only if not adding to existing
        plot(0, 0, type = "n",
             xlim = range(grid$x),
             ylim = range(grid$y),
             xlab = xlab, ylab = ylab,
             las = 1,
             axes = axes,
             main = main)
        if (!axes) {
            box()
        }
    }

    ## Plot maxima
    if (nrow(critical.points$maxima) > 0) {
        points(grid$x[critical.points$maxima[,1]],
               grid$y[critical.points$maxima[,2]],
               col = max.color, pch = max.pch, cex = point.size)
    }

    ## Plot minima
    if (nrow(critical.points$minima) > 0) {
        points(grid$x[critical.points$minima[,1]],
               grid$y[critical.points$minima[,2]],
               col = min.color, pch = min.pch, cex = point.size)
    }

    ## Return invisibly for potential additions to plot
    invisible(list(max.color = max.color,
                  min.color = min.color))
}

plot.gradient.trajectory <- function(grid, trajectory,
                                   ascending.color = "orange",
                                   descending.color = "purple",
                                   line.width = 2,
                                   add = FALSE) {
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
             main = "Gradient Trajectory")
    }

    ## Plot ascending path
    lines(asc.x, asc.y,
          col = ascending.color,
          lwd = line.width)

    ## Plot descending path
    lines(desc.x, desc.y,
          col = descending.color,
          lwd = line.width)
}

plot.function.contours <- function(grid, f.grid,
                                 main = "Function Contours",
                                 nlevels = 20,
                                 xlab = "", ylab = "",
                                 axes = FALSE,
                                 add = FALSE) {
    if (!add) {
        ## Create new contour plot
        contour(grid$x, grid$y, f.grid,
                nlevels = nlevels,
                xlab = xlab, ylab = ylab,
                las = 1,
                axes = axes,
                main = main)
    } else {
        ## Add contours to existing plot
        contour(grid$x, grid$y, f.grid,
                nlevels = nlevels,
                las = 1,
                axes = axes,
                add = TRUE)
    }
}

## Combined plotting function
plot.morse.analysis <- function(grid, f.grid, trajectory, critical.points,
                              show.contours = TRUE,
                              main = "Morse Analysis") {
    ## Start with contours if requested
    if (show.contours) {
        plot.function.contours(grid, f.grid, main = main)
    } else {
        ## If no contours, initialize empty plot
        plot(0, 0, type = "n",
             xlim = range(grid$x),
             ylim = range(grid$y),
             xlab = "x", ylab = "y",
             main = main)
    }

    ## Add critical points
    colors <- plot.grid.critical.points(grid, critical.points, add = TRUE)

    ## Add trajectory
    plot.gradient.trajectory(grid, trajectory,
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


find.ascending.trajectory.from.point <- function(x, y, step, mixture) {

    # Get gradient at current point
    grad <- mixture$gradient(x, y)

    # If gradient magnitude is very small, we're at/near a critical point
    if (sqrt(sum(grad^2)) < 1e-6) return(NULL)

    # Initialize variables for finding best neighbor
    max_alignment <- -Inf
    best_idx <- NULL

    # Check each neighbor
    for (k in 1:nrow(neighbors)) {
        ni <- neighbors[k,1]
        nj <- neighbors[k,2]
        # Convert neighbor indices to actual coordinates
        neighbor_point <- c(grid$x[ni], grid$y[nj])

        # Compute direction vector to neighbor
        direction <- neighbor_point - current_point
        direction_norm <- sqrt(sum(direction^2))

        if (direction_norm > 0) {
            # Compute alignment of direction with gradient (dot product)
            alignment <- sum(direction * grad) / direction_norm

            if (alignment > max_alignment) {
                max_alignment <- alignment
                best_idx <- k
            }
        }
    }

    # Return NULL if no good ascending direction found
    if (max_alignment <= 0 || is.null(best_idx)) return(NULL)
    return(neighbors[best_idx,])
}


find.next.ascending.point <- function(i, j, grid, mixture) {
    neighbors <- get.neighbors(i, j, length(grid$x))
    current_point <- c(grid$x[i], grid$y[j])

    # Get gradient at current point
    grad <- mixture$gradient(current_point[1], current_point[2])

    # If gradient magnitude is very small, we're at/near a critical point
    if (sqrt(sum(grad^2)) < 1e-6) return(NULL)

    # Initialize variables for finding best neighbor
    max_alignment <- -Inf
    best_idx <- NULL

    # Check each neighbor
    for (k in 1:nrow(neighbors)) {
        ni <- neighbors[k,1]
        nj <- neighbors[k,2]
        # Convert neighbor indices to actual coordinates
        neighbor_point <- c(grid$x[ni], grid$y[nj])

        # Compute direction vector to neighbor
        direction <- neighbor_point - current_point
        direction_norm <- sqrt(sum(direction^2))

        if (direction_norm > 0) {
            # Compute alignment of direction with gradient (dot product)
            alignment <- sum(direction * grad) / direction_norm

            if (alignment > max_alignment) {
                max_alignment <- alignment
                best_idx <- k
            }
        }
    }

    # Return NULL if no good ascending direction found
    if (max_alignment <= 0 || is.null(best_idx)) return(NULL)
    return(neighbors[best_idx,])
}

find.next.descending.point <- function(i, j, grid, mixture) {
    neighbors <- get.neighbors(i, j, length(grid$x))
    current_point <- c(grid$x[i], grid$y[j])

    # Get gradient at current point
    grad <- mixture$gradient(current_point[1], current_point[2])

    # If gradient magnitude is very small, we're at/near a critical point
    if (sqrt(sum(grad^2)) < 1e-6) return(NULL)

    # Initialize variables for finding best neighbor
    max_alignment <- -Inf
    best_idx <- NULL

    # Check each neighbor
    for (k in 1:nrow(neighbors)) {
        ni <- neighbors[k,1]
        nj <- neighbors[k,2]
        # Convert neighbor indices to actual coordinates
        neighbor_point <- c(grid$x[ni], grid$y[nj])

        # Compute direction vector to neighbor
        direction <- neighbor_point - current_point
        direction_norm <- sqrt(sum(direction^2))

        if (direction_norm > 0) {
            # Compute alignment with negative gradient (for descent)
            alignment <- sum(direction * (-grad)) / direction_norm

            if (alignment > max_alignment) {
                max_alignment <- alignment
                best_idx <- k
            }
        }
    }

    # Return NULL if no good descending direction found
    if (max_alignment <= 0 || is.null(best_idx)) return(NULL)
    return(neighbors[best_idx,])
}

compute.gradient.trajectory <- function(i, j, grid, mixture) {
    # Initialize paths
    ascending.path <- matrix(c(i, j), ncol = 2)
    descending.path <- matrix(c(i, j), ncol = 2)

    # Compute ascending path
    current <- c(i, j)
    while (!is.null(current)) {
        next.point <- find.next.ascending.point(current[1], current[2],
                                                grid, mixture)
        if (!is.null(next.point)) {
            ascending.path <- rbind(ascending.path, next.point)
            current <- next.point
        } else {
            break
        }
    }
    local.max <- current

    # Compute descending path
    current <- c(i, j)
    while (!is.null(current)) {
        next.point <- find.next.descending.point(current[1], current[2],
                                                 grid, mixture)
        if (!is.null(next.point)) {
            descending.path <- rbind(descending.path, next.point)
            current <- next.point
        } else {
            break
        }
    }
    local.min <- current

    return(list(
        ascending.path = ascending.path,
        descending.path = descending.path,
        local.max = local.max,
        local.min = local.min
    ))
}

compute.all.trajectories <- function(grid, mixture) {
    n <- length(grid$x)
    points.to.process <- expand.grid(i = 1:n, j = 1:n)
    trajectories <- list()

    while(nrow(points.to.process) > 0) {
        cat("\rnrow(points.to.process): ", nrow(points.to.process))
        # Take first point
        current <- points.to.process[1,]

        # Compute its trajectory
        traj <- compute.gradient.trajectory(current$i, current$j,
                                            grid, mixture)

        # Remove all points in the trajectory from points.to.process
        all.points <- rbind(traj$ascending.path, traj$descending.path)
        points.to.process <- points.to.process[!apply(points.to.process, 1,
            function(p) any(apply(all.points, 1,
                function(q) all(p == q)))), ]

        # Add trajectory to list
        trajectories[[length(trajectories) + 1]] <- traj
    }

    return(trajectories)
}

compute.morse.smale.cells <- function(trajectories) {
    # First, let's identify all unique local maxima and minima
    local_maxima <- unique(do.call(rbind, lapply(trajectories, function(t) t$local.max)))
    local_minima <- unique(do.call(rbind, lapply(trajectories, function(t) t$local.min)))

    # Create a lookup table for cell assignments
    # Each cell will be identified by a pair of indices (max_idx, min_idx)
    cell_lookup <- list()

    # For each trajectory, create a cell identifier and assign all points
    # in the trajectory to that cell
    for (traj_idx in 1:length(trajectories)) {
        traj <- trajectories[[traj_idx]]

        # Find indices of this trajectory's extrema in our unique lists
        max_idx <- which(apply(local_maxima, 1, function(x)
            all(x == traj$local.max)))
        min_idx <- which(apply(local_minima, 1, function(x)
            all(x == traj$local.min)))

        # Create cell identifier
        cell_id <- paste(max_idx, min_idx, sep = "_")

        # Get all points in this trajectory (excluding the extrema themselves)
        traj_points <- rbind(
            traj$ascending.path[-nrow(traj$ascending.path),],  # exclude maximum
            traj$descending.path[-nrow(traj$descending.path),] # exclude minimum
        )

        # Add these points to the cell lookup
        if (is.null(cell_lookup[[cell_id]])) {
            cell_lookup[[cell_id]] <- traj_points
        } else {
            cell_lookup[[cell_id]] <- unique(rbind(cell_lookup[[cell_id]],
                                                 traj_points))
        }
    }

    # Create a more structured output
    morse_smale_complex <- list(
        cells = cell_lookup,
        local_maxima = local_maxima,
        local_minima = local_minima,
        cell_summary = data.frame(
            cell_id = names(cell_lookup),
            max_idx = sapply(strsplit(names(cell_lookup), "_"), `[`, 1),
            min_idx = sapply(strsplit(names(cell_lookup), "_"), `[`, 2),
            point_count = sapply(cell_lookup, nrow)
        )
    )

    # Add helper function to find which cell contains a point
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



plot.morse.smale.cells <- function(grid, morse_smale_complex, f.grid = NULL) {
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
    cell_colors <- hcl.colors(n.cells, "Set 2")  # Using HCL colors for better distinction

    # Create a named vector for all possible labels
    color_mapping <- c(
        setNames(cell_colors, unique_cells),
        local_maximum = "red",
        local_minimum = "blue"
    )

    # Set up the plotting area with appropriate margins
    par(mar = c(5, 4, 4, 6))  # Adjust margins to accommodate legend

    # Create base plot
    image(grid$x, grid$y, cell.labels,
          col = color_mapping,
          xlab = "x",
          ylab = "y",
          main = "Morse-Smale Complex")

    # Add contour lines if f.grid is provided
    if (!is.null(f.grid)) {
        contour(grid$x, grid$y, f.grid,
                add = TRUE,
                drawlabels = FALSE,
                col = "black",
                lty = 2)
    }

    # Add cell boundaries
    contour(grid$x, grid$y, cell.labels,
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
           fill = legend_colors,
           pch = legend_pch,
           pt.bg = legend_colors,
           title = "Morse-Smale Components",
           inset = c(-0.2, 0),
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
        grad[1] <- 0  # Zero out x component
    }
    if(at_y_min || at_y_max) {
        grad[2] <- 0  # Zero out y component
    }

    return(grad)
}

find.ascending.trajectory.from.point <- function(x, y, step = 0.01, mixture,
                                                 max_steps = 1000,
                                                 grad_threshold = 1e-5,
                                                 angle_threshold = -0.5,
                                                 domain = list(x_min = 0, x_max = 1,
                                                               y_min = 0, y_max = 1)) {
    # Initialize trajectory storage
    trajectory <- matrix(nrow = max_steps, ncol = 2)
    trajectory[1,] <- c(x, y)

    # Keep track of previous gradient for angle comparison
    prev_grad <- mixture$gradient(x, y)
    prev_grad_norm <- sqrt(sum(prev_grad^2))
    if(prev_grad_norm > 0) {
        prev_grad <- prev_grad / prev_grad_norm
    }

    n_points <- 1

    for(i in 2:max_steps) {
        # Get gradient at current point
        grad <- mixture$gradient(x, y)

        # Project gradient if we're on boundary
        grad <- project.gradient.to.boundary(x, y, grad, domain)

        grad_norm <- sqrt(sum(grad^2))

        # Stop if gradient is too small (we're at a critical point or stuck at corner)
        if(grad_norm < grad_threshold) {
            break
        }

        # Normalize gradient
        grad <- grad / grad_norm

        # Check angle between current and previous gradient
        angle_cos <- sum(grad * prev_grad)
        if(angle_cos < angle_threshold) {
            break
        }

        # Update position
        new_x <- x + step * grad[1]
        new_y <- y + step * grad[2]

        # Ensure we stay within domain
        x <- min(max(new_x, domain$x_min), domain$x_max)
        y <- min(max(new_y, domain$y_min), domain$y_max)

        # Store new position
        trajectory[i,] <- c(x, y)
        n_points <- i

        # Update previous gradient
        prev_grad <- grad
    }

    return(trajectory[1:n_points,])
}


## find.ascending.trajectory.from.point <- function(x, y, step = 0.01, mixture, max_steps = 1000,
##                                                grad_threshold = 1e-5, angle_threshold = -0.5) {
##     # Initialize trajectory storage
##     trajectory <- matrix(nrow = max_steps, ncol = 2)
##     trajectory[1,] <- c(x, y)

##     # Keep track of previous gradient for angle comparison
##     prev_grad <- mixture$gradient(x, y)
##     prev_grad_norm <- sqrt(sum(prev_grad^2))
##     if(prev_grad_norm > 0) {
##         prev_grad <- prev_grad / prev_grad_norm
##     }

##     n_points <- 1

##     for(i in 2:max_steps) {
##         # Get gradient at current point
##         grad <- mixture$gradient(x, y)
##         grad_norm <- sqrt(sum(grad^2))

##         # Stop if gradient is too small (we're at a critical point)
##         if(grad_norm < grad_threshold) {
##             break
##         }

##         # Normalize gradient
##         grad <- grad / grad_norm

##         # Check angle between current and previous gradient
##         # cos(angle) = dot product of normalized vectors
##         angle_cos <- sum(grad * prev_grad)
##         if(angle_cos < angle_threshold) {
##             # Gradient direction changed significantly, likely crossed a ridge
##             break
##         }

##         # Update position
##         x <- x + step * grad[1]
##         y <- y + step * grad[2]

##         # Store new position
##         trajectory[i,] <- c(x, y)
##         n_points <- i

##         # Update previous gradient
##         prev_grad <- grad
##     }

##     # Return only the valid part of trajectory
##     return(trajectory[1:n_points,])
## }

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

# Function to plot trajectories and critical points
plot.morse.smale.complex <- function(grid, mixture, n_sample_points = 20) {
    # Create base contour plot
    f_grid <- evaluate.function.on.grid(mixture$f, grid)
    plot.function.contours(grid, f_grid, main = "Morse-Smale Complex", nlevels = 10)

    # Plot trajectories from sample points
    x_vals <- seq(grid$x_min, grid$x_max, length.out = n_sample_points)
    y_vals <- seq(grid$y_min, grid$y_max, length.out = n_sample_points)

    for(x in x_vals) {
        for(y in y_vals) {
            # Find ascending and descending trajectories
            asc_traj <- find.ascending.trajectory.from.point(x, y, 0.01, mixture)
            desc_traj <- find.descending.trajectory.from.point(x, y, 0.01, mixture)

            # Plot trajectories
            lines(asc_traj[,1], asc_traj[,2], col = "red", lwd = 0.5)
            lines(desc_traj[,1], desc_traj[,2], col = "blue", lwd = 0.5)
        }
    }

    # Add critical points
    plot.grid.critical.points(grid, mixture$critical.points, add = TRUE)
}

# Helper function to generate points on a circle, restricted to domain if needed
generate.circle.points <- function(center_x, center_y, radius, n_points,
                                   domain = list(x_min = 0, x_max = 1,
                                                 y_min = 0, y_max = 1)) {
    ## Generate angles uniformly around the circle
    large.n <- 10000
    angles <- seq(0, 2*pi, length.out = large.n + 1)[-(large.n + 1)]

    # Generate points on the circle
    points <- matrix(nrow = large.n, ncol = 2)
    points[,1] <- center_x + radius * cos(angles)
    points[,2] <- center_y + radius * sin(angles)

    # Filter points to keep only those within domain
    valid_points <- points[,1] >= domain$x_min & points[,1] <= domain$x_max & points[,2] >= domain$y_min & points[,2] <= domain$y_max
    angles <- angles[valid_points]
    angles <- seq(min(angles), max(angles), length.out = n_points + 1)[-(n_points + 1)]

    points <- matrix(nrow = n_points, ncol = 2)
    points[,1] <- center_x + radius * cos(angles)
    points[,2] <- center_y + radius * sin(angles)

    return(points)
}

plot.morse.smale.complex.from.critical <- function(grid, mixture,
                                                 critical.points,
                                                 epsilon = 0.05,
                                                 n_circle_points = 32) {
    # Create base contour plot
    f_grid <- evaluate.function.on.grid(mixture$f, grid)
    plot.function.contours(grid, f_grid, main = "Morse-Smale Complex", nlevels = 10)

    # For each critical point
    for(i in 1:nrow(critical.points)) {
        # Generate circle points around the critical point
        circle_points <- generate.circle.points(
            critical.points$x[i],
            critical.points$y[i],
            epsilon,
            n_circle_points,
            grid
        )

        # If this is a maximum, compute descending trajectories
        # If this is a minimum, compute ascending trajectories
        if(critical.points$type[i] == "maximum") {
            for(j in 1:nrow(circle_points)) {
                traj <- find.descending.trajectory.from.point(
                    circle_points[j,1],
                    circle_points[j,2],
                    0.01,
                    mixture
                )
                lines(traj[,1], traj[,2], col = "blue", lwd = 0.5)
            }
        } else {  # minimum
            for(j in 1:nrow(circle_points)) {
                traj <- find.ascending.trajectory.from.point(
                    circle_points[j,1],
                    circle_points[j,2],
                    0.01,
                    mixture
                )
                lines(traj[,1], traj[,2], col = "red", lwd = 0.5)
            }
        }
    }

    # Add critical points
    plot.grid.critical.points(grid, critical.points, add = TRUE)
}

plot.morse.smale.trajectories.from.point <- function(grid, mixture,
                                                    point_x, point_y,
                                                    is_maximum = TRUE,
                                                    epsilon = 0.05,
                                                    n_circle_points = 32) {
    # Generate circle points
    circle_points <- generate.circle.points(point_x, point_y, epsilon,
                                            n_circle_points, grid)

    # Plot base contours if this is the first plot
    if(!par("new")) {
        f_grid <- evaluate.function.on.grid(mixture$f, grid)
        plot.function.contours(grid, f_grid,
                             main = "Gradient Flow from Critical Point",
                             nlevels = 10)
    }

    # Compute and plot trajectories
    for(i in 1:nrow(circle_points)) {
        if(is_maximum) {
            traj <- find.descending.trajectory.from.point(
                circle_points[i,1],
                circle_points[i,2],
                0.01,
                mixture
            )
            lines(traj[,1], traj[,2], col = "blue", lwd = 0.5)
        } else {
            traj <- find.ascending.trajectory.from.point(
                circle_points[i,1],
                circle_points[i,2],
                0.01,
                mixture
            )
            lines(traj[,1], traj[,2], col = "red", lwd = 0.5)
        }
    }

    # Add the critical point
    points(point_x, point_y, pch = 19,
           col = if(is_maximum) "red" else "blue", cex = 1.5)
}



find.critical.points.improved <- function(grid, mixture, n_trajectory_points = 20) {

    f.grid <- evaluate.function.on.grid(mixture, grid)

    ## Step 1: Grid-based initial detection
    initial_points <- find.grid.critical.points(grid, f.grid)

    # Step 2: Refinement using trajectories
    refined_points <- refine.critical.points(initial_points, mixture, n_trajectory_points)

    return(refined_points)
}

## find.grid.critical.points <- function(grid, f) {
##     # Evaluate function on grid
##     z <- matrix(0, length(grid$x), length(grid$y))
##     for(i in 1:length(grid$x)) {
##         for(j in 1:length(grid$y)) {
##             z[i,j] <- f(grid$x[i], grid$y[j])
##         }
##     }

##     # Initialize lists for different types of critical points
##     maxima <- matrix(ncol = 2, nrow = 0)
##     minima <- matrix(ncol = 2, nrow = 0)
##     saddles <- matrix(ncol = 2, nrow = 0)

##     # Check interior points
##     for(i in 2:(length(grid$x)-1)) {
##         for(j in 2:(length(grid$y)-1)) {
##             neighborhood <- z[i-1:i+1, j-1:j+1]
##             center <- z[i,j]

##             if(all(center >= neighborhood) && any(center > neighborhood)) {
##                 maxima <- rbind(maxima, c(grid$x[i], grid$y[j]))
##             } else if(all(center <= neighborhood) && any(center < neighborhood)) {
##                 minima <- rbind(minima, c(grid$x[i], grid$y[j]))
##             }
##         }
##     }

##     # Check boundaries
##     check_boundary_points(grid, z, maxima, minima)

##     return(list(maxima = maxima, minima = minima))
## }

refine.critical.points <- function(initial_points, mixture, n_points = 10, radius = 0.02) {
    refined_maxima <- matrix(ncol = 2, nrow = 0)

    # For each maximum, generate circle of points and follow trajectories
    for(i in 1:nrow(initial_points$coordinates$maxima)) {
        x <- initial_points$coordinates$maxima[i,1]
        y <- initial_points$coordinates$maxima[i,2]

        # Generate circle of points around the maximum
        circle_points <- generate.circle.points(x, y, radius = radius, n_points)

        # Follow trajectories and collect endpoints
        endpoints <- matrix(ncol = 2, nrow = n_points)
        for(j in 1:n_points) {
            trajectory <- find.ascending.trajectory.from.point(
                circle_points[j,1], circle_points[j,2],
                step = 0.01, mixture
            )
            endpoints[j,] <- trajectory[nrow(trajectory),]
        }

        # Compute refined position as median of endpoints
        refined_pos <- apply(endpoints, 2, median)
        refined_maxima <- rbind(refined_maxima, refined_pos)

        # Add uncertainty estimate
        attr(refined_maxima, "uncertainty") <- apply(endpoints, 2, sd)
    }

    return(list(
        maxima = refined_maxima,
        initial_points = initial_points
    ))
}

#' Add Arrow to Trajectory Line with Controllable Orientation
#'
#' This function draws a trajectory line with an arrow that can be oriented either
#' automatically (following the trajectory) or manually (using a specified angle).
#'
#' @param trajectory Matrix (n x 2) containing trajectory coordinates
#' @param p Numeric between 0 and 1 indicating relative position of arrow (default: 0.5)
#' @param arrow.length Numeric specifying arrow head length (default: 0.02)
#' @param arrow.angle Numeric specifying arrow head angle in degrees (default: 30)
#' @param line.col Color for the trajectory line (default: "gray")
#' @param line.lwd Line width
#' @param arrow.col Color for the arrow head (default: same as line.col)
#' @param orientation Either "auto" to follow trajectory direction, or a number
#'        specifying the angle in degrees (0 = right, 90 = up, etc.)
#' @param reverse Logical indicating whether to reverse the arrow direction (default: FALSE)
plot.trajectory.with.arrow <- function(trajectory,
                                     p = 0.5,
                                     arrow.length = 0.02,
                                     arrow.angle = 30,
                                     line.col = "gray",
                                     line.lwd = 2,
                                     arrow.col = line.col,
                                     orientation = "auto",
                                     reverse = FALSE) {
    # Draw the trajectory line
    lines(trajectory, col = line.col, lwd = line.lwd)

    # Calculate the position along the trajectory
    n_points <- nrow(trajectory)
    distances <- sqrt(diff(trajectory[,1])^2 + diff(trajectory[,2])^2)
    cum_distances <- c(0, cumsum(distances))
    total_distance <- sum(distances)
    target_distance <- p * total_distance

    # Find the segment containing our target position
    segment_idx <- max(which(cum_distances <= target_distance))

    if (segment_idx < n_points) {
        # Calculate exact arrow position through interpolation
        remaining_distance <- target_distance - cum_distances[segment_idx]
        segment_length <- distances[segment_idx]
        t <- remaining_distance / segment_length

        arrow_pos <- c(
            trajectory[segment_idx,1] + t * (trajectory[segment_idx + 1,1] - trajectory[segment_idx,1]),
            trajectory[segment_idx,2] + t * (trajectory[segment_idx + 1,2] - trajectory[segment_idx,2])
        )

        # Determine arrow orientation
        if (orientation == "auto") {
            # Use trajectory direction
            dx <- trajectory[segment_idx + 1,1] - trajectory[segment_idx,1]
            dy <- trajectory[segment_idx + 1,2] - trajectory[segment_idx,2]
            angle <- atan2(dy, dx)
        } else {
            # Use specified orientation (convert degrees to radians)
            angle <- orientation * pi / 180
        }

        # Reverse direction if requested
        if (reverse) {
            angle <- angle + pi
        }

        # Calculate arrow head points
        arrow_angle_rad <- arrow.angle * pi / 180
        left_angle <- angle + arrow_angle_rad
        right_angle <- angle - arrow_angle_rad

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
}

identify.destination.minimum <- function(trajectory, minima_coordinates, tolerance = 0.05) {
    # Get the final point of the trajectory
    final_point <- trajectory[nrow(trajectory), ]

    # Calculate distances from the final point to all minima
    distances <- apply(minima_coordinates, 1, function(min_coord) {
        sqrt(sum((final_point - min_coord)^2))
    })

    # If the smallest distance is larger than our tolerance, warn about potential issues
    min_distance <- min(distances)
    if (min_distance > tolerance) {
        warning(sprintf("Trajectory endpoint (%.3f, %.3f) is relatively far (%.3f) from nearest minimum. Consider adjusting tolerance.",
                       final_point[1], final_point[2], min_distance))
    }

    # Find which minimum is closest
    closest_minimum_index <- which.min(distances)

    # Return the label for this minimum
    return(paste0("m", closest_minimum_index))
}


# Function to find maximum along a domain edge
find.edge.maximum <- function(f1, edge, n.points = 100) {
    # edge should be one of: "bottom", "top", "left", "right"
    # Returns the position and value of maximum along specified edge

    if (edge %in% c("bottom", "top")) {
        # For horizontal edges, vary x from 0 to 1
        x.vals <- seq(0, 1, length.out = n.points)
        y.val <- ifelse(edge == "bottom", 0, 1)

        # Evaluate function along the edge
        vals <- sapply(x.vals, function(x) f1$f(x, y.val))
        max.idx <- which.max(vals)

        return(list(
            position = c(x.vals[max.idx], y.val),
            value = vals[max.idx]
        ))
    } else {
        # For vertical edges, vary y from 0 to 1
        y.vals <- seq(0, 1, length.out = n.points)
        x.val <- ifelse(edge == "left", 0, 1)

        # Evaluate function along the edge
        vals <- sapply(y.vals, function(y) f1$f(x.val, y))
        max.idx <- which.max(vals)

        return(list(
            position = c(x.val, y.vals[max.idx]),
            value = vals[max.idx]
        ))
    }
}

# Function to find all boundary maxima
find.all.boundary.maxima <- function(f1) {
    edges <- c("bottom", "top", "left", "right")
    maxima <- lapply(edges, function(edge) find.edge.maximum(f1, edge))
    names(maxima) <- edges
    return(maxima)
}

find.interior.saddle <- function(f1, M1.pos, M2.pos, n.points = 100) {
    # Create points along the line from M1 to M2
    t <- seq(0, 1, length.out = n.points)
    points <- t(sapply(t, function(s) {
        M1.pos * (1-s) + M2.pos * s  # Linear interpolation
    }))

    # Evaluate function along this line
    vals <- apply(points, 1, function(x) f1$f(x[1], x[2]))

    # Find the minimum (this will be our saddle)
    min.idx <- which.min(vals)

    return(list(
        position = points[min.idx,],
        value = vals[min.idx]
    ))
}

compute.separatrices <- function(f1, boundary_maxima, interior_saddle,
                               step = 0.01, max_steps = 1000) {
    # This function computes the separating curves (separatrices) that divide our domain
    # into Morse-Smale cells. It uses the exact gradient information from our Gaussian mixture.

    separatrices <- list()

    # Helper function to compute a single trajectory
    # This function follows the gradient flow either up (ascending) or down (descending)
    compute_trajectory <- function(start_point, ascending = FALSE) {
        trajectory <- matrix(NA, nrow = max_steps, ncol = 2)
        trajectory[1,] <- start_point

        for (i in 2:max_steps) {
            # Get the exact gradient at our current position
            grad <- f1$gradient(trajectory[i-1, 1], trajectory[i-1, 2])

            # For ascending trajectories, we follow the gradient
            # For descending trajectories, we follow the negative gradient
            if (!ascending) grad <- -grad

            # Normalize the gradient to ensure consistent step sizes
            grad_norm <- sqrt(sum(grad^2))
            if (grad_norm < 1e-10) break  # Stop if gradient is essentially zero

            grad <- grad / grad_norm
            new_point <- trajectory[i-1,] + step * grad

            # Stop if we've left the domain
            if (any(new_point < 0) || any(new_point > 1)) break

            trajectory[i,] <- new_point
        }

        # Remove unused rows and return only the valid part of the trajectory
        trajectory <- trajectory[!is.na(trajectory[,1]),]
        return(trajectory)
    }

    # Compute descending trajectories from each boundary maximum
    # These form part of the boundaries between Morse-Smale cells
    for (edge in names(boundary_maxima)) {
        start_point <- boundary_maxima[[edge]]$position
        separatrices[[paste0("boundary_", edge)]] <-
            compute_trajectory(start_point, ascending = FALSE)
    }

    # Compute trajectories from the interior saddle point
    # These are particularly important as they separate the domains of attraction
    saddle_point <- interior_saddle$position

    # For the saddle point, we need both ascending and descending trajectories
    # We'll use small perturbations in the eigendirections of the Hessian
    # but for now we'll use simple offsets
    eps <- step/10  # Small perturbation

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

plot.separatrices <- function(separatrices, M1_pos, M2_pos,
                              line.lwd = 2,
                              line.lty = 2,
                              line.col = "black",
                              add = TRUE) {
    ## Create new plot if needed
    if (!add) {
        plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             xlab = "", ylab = "")
    }

    ## For boundary separatrices, we connect them to the appropriate maximum
    ## The pattern we observe is:
    ## - Points on bottom and left edges connect to M1
    ## - Points on top and right edges connect to M2

    ## Connect bottom and left boundary maxima to M1
    for (edge in c("boundary_bottom", "boundary_left")) {
        if (!is.null(separatrices[[edge]])) {
            ## Draw line segment from boundary maximum to M1
            segments(x0 = separatrices[[edge]][1],  # x-coordinate of boundary point
                     y0 = separatrices[[edge]][2],  # y-coordinate of boundary point
                     x1 = M1_pos[1],                # x-coordinate of M1
                     y1 = M1_pos[2],                # y-coordinate of M1,
                     lty = line.lty,
                     lwd = line.lwd,
                     col = line.col)
        }
    }

    ## Connect top and right boundary maxima to M2
    for (edge in c("boundary_top", "boundary_right")) {
        if (!is.null(separatrices[[edge]])) {
            ## Draw line segment from boundary maximum to M2
            segments(x0 = separatrices[[edge]][1],  # x-coordinate of boundary point
                     y0 = separatrices[[edge]][2],  # y-coordinate of boundary point
                     x1 = M2_pos[1],                # x-coordinate of M2
                     y1 = M2_pos[2],                # y-coordinate of M2,
                     lty = line.lty,
                     lwd = line.lwd,
                     col = line.col)
        }
    }

    ## Plot interior saddle separatrices as before
    ## for (traj_type in c("saddle_ascending1", "saddle_ascending2",
    ##                    "saddle_descending1", "saddle_descending2")) {
    ##     trajectory <- separatrices[[traj_type]]
    ##     if (!is.null(trajectory) && is.matrix(trajectory)) {
    ##         lines(trajectory[, 1], trajectory[, 2],
    ##               lty = line.lty,
    ##               lwd = line.lwd,
    ##               col = line.col)
    ##     }
    ## }

    ## Draw a line connecting M1 and M2
    segments(x0 = M1_pos[1],                # x-coordinate of M1
             y0 = M1_pos[2],                # y-coordinate of M1,
             x1 = M2_pos[1],                # x-coordinate of M2
             y1 = M2_pos[2],                # y-coordinate of M2,
             lty = line.lty,
             lwd = line.lwd,
             col = line.col)

    ## Draw a line connecting the end of separatrices$saddle_descending1 with separatrices$saddle_descending2
    d1 <- separatrices$saddle_descending1
    d1.pt <- d1[nrow(d1),]
    d2 <- separatrices$saddle_descending2
    d2.pt <- d2[nrow(d2),]
    segments(x0 = d1.pt[1],                 # x-coordinate of d1.pt
             y0 = d1.pt[2],                 # y-coordinate of d1.pt
             x1 = d2.pt[1],                 # x-coordinate of d2.pt
             y1 = d2.pt[2],                 # y-coordinate of d2.pt
             lty = line.lty,
             lwd = line.lwd,
             col = line.col)

}


plot.separatrices.with.cells <- function(separatrices, M1.pos, M2.pos, saddle_pos, add = TRUE, cell_colors = list(
        "M1_m1" = rgb(0, 0.5, 0, alpha = 0.1),    # Light green
        "M1_m2" = rgb(1, 1, 0, alpha = 0.1),      # Light yellow
        "M1_m3" = rgb(1, 0, 0, alpha = 0.1),      # Light red
        "M2_m2" = rgb(0, 0, 1, alpha = 0.1),      # Light blue
        "M2_m3" = rgb(0.6, 0, 1, alpha = 0.1),    # Light purple
        "M2_m4" = rgb(0, 1, 1, alpha = 0.1)       # Light turquoise
        )) {

    # Create new plot if needed
    if (!add) {
        plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             xlab = "", ylab = "")
    }

    ## Fill cells using polygons

    ## M1_m1 cell
    polygon(
        x = c(0, separatrices$boundary_left[1], M1.pos[1], separatrices$boundary_bottom[1], 0),
        y = c(0, separatrices$boundary_left[2], M1.pos[2], separatrices$boundary_bottom[2], 0),
        col = cell_colors$M1_m1,
        border = NA  # No border for the filled region
    )

    ## Finding the saddle point as the interesection of d1.pt-d2.pt and M1.pos-M2.pos line segments
    d1 <- separatrices$saddle_descending1
    d1.pt <- d1[nrow(d1),]
    d2 <- separatrices$saddle_descending2
    d2.pt <- d2[nrow(d2),]
    saddle.pt <- find.line.intersection(d1.pt, d2.pt, M1.pos, M2.pos)

    ## M1_m2 cell - technically it also contains a line segment [d2.pt, m2], but we ignore it
    ## d2.pt -> saddle.pt -> M1 -> separatrices$boundary_left -> d2.pt
    polygon(
        x = c(d2.pt[1], saddle.pt[1], M1.pos[1], separatrices$boundary_left[1], d2.pt[1]),
        y = c(d2.pt[2], saddle.pt[2], M1.pos[2], separatrices$boundary_left[2], d2.pt[2]),
        col = cell_colors$M1_m2,
        border = NA  # No border for the filled region
    )

    ## M1_m3 cell - technically it also contains a line segment [d2.pt, m3], but we ignore it
    ## d1.pt -> saddle.pt -> M1 -> separatrices$boundary_bottom -> d1.pt
    polygon(
        x = c(d1.pt[1], saddle.pt[1], M1.pos[1], separatrices$boundary_bottom[1], d1.pt[1]),
        y = c(d1.pt[2], saddle.pt[2], M1.pos[2], separatrices$boundary_bottom[2], d1.pt[2]),
        col = cell_colors$M1_m3,
        border = NA  # No border for the filled region
    )

    ## M2_m4 cell
    polygon(
        x = c(1, separatrices$boundary_top[1], M2.pos[1], separatrices$boundary_right[1], 1),
        y = c(1, separatrices$boundary_top[2], M2.pos[2], separatrices$boundary_right[2], 1),
        col = cell_colors$M2_m4,
        border = NA  # No border for the filled region
    )

    ## M2_m2 cell - technically it also contains a line segment [d2.pt, m2], but we ignore it
    ## m2 -> separatrices$boundary_top -> M2 -> saddle.pt -> d2.pt -> m2
    polygon(
        x = c(0, separatrices$boundary_top[1], M2.pos[1], saddle.pt[1], d2.pt[1], 0),
        y = c(1, separatrices$boundary_top[2], M2.pos[2], saddle.pt[2], d2.pt[2], 1),
        col = cell_colors$M2_m2,
        border = NA  # No border for the filled region
    )

    ## M2_m3 cell - technically it also contains a line segment [d2.pt, m3], but we ignore it
    ## m3 -> d1.pt -> saddle.pt -> M2 -> separatrices$boundary_right -> m3
    polygon(
        x = c(1, d1.pt[1], saddle.pt[1], M2.pos[1], separatrices$boundary_right[1], 1),
        y = c(0, d1.pt[2], saddle.pt[2], M2.pos[2], separatrices$boundary_right[2], 0),
        col = cell_colors$M2_m3,
        border = NA  # No border for the filled region
    )
}


v2.plot.separatrices.with.cells <- function(separatrices, M1.pos, M2.pos, saddle_pos, add = TRUE, cell_colors = list(
        "M1_m1" = rgb(0, 0.5, 0, alpha = 0.1),    # Light green
        "M1_m2" = rgb(1, 1, 0, alpha = 0.1),      # Light yellow
        "M1_m3" = rgb(1, 0, 0, alpha = 0.1),      # Light red
        "M2_m2" = rgb(0, 0, 1, alpha = 0.1),      # Light blue
        "M2_m3" = rgb(0.6, 0, 1, alpha = 0.1),    # Light purple
        "M2_m4" = rgb(0, 1, 1, alpha = 0.1)       # Light turquoise
        )) {

    # Create new plot if needed
    if (!add) {
        plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             xlab = "", ylab = "")
    }

    ## Fill cells using polygons

    ## M1_m1 cell
    polygon(
        x = c(0, separatrices$boundary_left[1], M1.pos[1], separatrices$boundary_bottom[1], 0),
        y = c(0, separatrices$boundary_left[2], M1.pos[2], separatrices$boundary_bottom[2], 0),
        col = cell_colors$M1_m1,
        border = NA  # No border for the filled region
    )

    ## Finding the saddle point as the interesection of d1.pt-d2.pt and M1.pos-M2.pos line segments
    d1 <- separatrices$saddle_descending1
    d1.pt <- d1[nrow(d1),]
    d2 <- separatrices$saddle_descending2
    d2.pt <- d2[nrow(d2),]
    saddle.pt <- find.line.intersection(d1.pt, d2.pt, M1.pos, M2.pos)

    ## M1_m2 cell - technically it also contains a line segment [d2.pt, m2], but we ignore it
    ## m2 -> d2.pt -> saddle.pt -> M1 -> separatrices$boundary_left -> d2.pt -> m2
    polygon(
        x = c(0, d2.pt[1], saddle.pt[1], M1.pos[1], separatrices$boundary_left[1], 0),
        y = c(1, d2.pt[2], saddle.pt[2], M1.pos[2], separatrices$boundary_left[2], 1),
        col = cell_colors$M1_m2,
        border = NA  # No border for the filled region
    )

    ## M1_m3 cell - technically it also contains a line segment [d2.pt, m3], but we ignore it
    ## d1.pt -> saddle.pt -> M1 -> separatrices$boundary_bottom -> m3 -> d1.pt
    polygon(
        x = c(d1.pt[1], saddle.pt[1], M1.pos[1], separatrices$boundary_bottom[1], 1, d1.pt[1]),
        y = c(d1.pt[2], saddle.pt[2], M1.pos[2], separatrices$boundary_bottom[2], 0, d1.pt[2]),
        col = cell_colors$M1_m3,
        border = NA  # No border for the filled region
    )

    ## M2_m4 cell
    polygon(
        x = c(1, separatrices$boundary_top[1], M2.pos[1], separatrices$boundary_right[1], 1),
        y = c(1, separatrices$boundary_top[2], M2.pos[2], separatrices$boundary_right[2], 1),
        col = cell_colors$M2_m4,
        border = NA  # No border for the filled region
    )

    ## M2_m2 cell - technically it also contains a line segment [d2.pt, m2], but we ignore it
    ## m2 -> separatrices$boundary_top -> M2 -> saddle.pt -> d2.pt -> m2
    polygon(
        x = c(0, separatrices$boundary_top[1], M2.pos[1], saddle.pt[1], d2.pt[1], 0),
        y = c(1, separatrices$boundary_top[2], M2.pos[2], saddle.pt[2], d2.pt[2], 1),
        col = cell_colors$M2_m2,
        border = NA  # No border for the filled region
    )

    ## M2_m3 cell - technically it also contains a line segment [d2.pt, m3], but we ignore it
    ## m3 -> d1.pt -> saddle.pt -> M2 -> separatrices$boundary_right -> m3
    polygon(
        x = c(1, d1.pt[1], saddle.pt[1], M2.pos[1], separatrices$boundary_right[1], 1),
        y = c(0, d1.pt[2], saddle.pt[2], M2.pos[2], separatrices$boundary_right[2], 0),
        col = cell_colors$M2_m3,
        border = NA  # No border for the filled region
    )
}



# Create a Point class using R's S3 system
create_point <- function(x, y) {
  point <- list(x = x, y = y)
  class(point) <- "Point"
  return(point)
}

#' Find the intersection point of two line segments
#'
#' This function calculates the intersection point of two line segments if it exists.
#' It uses the parametric form of line segments and solves the system of equations:
#' P1 + t(P2-P1) = Q1 + s(Q2-Q1)
#' where P1,P2 are points of first line segment and Q1,Q2 are points of second line segment.
#'
#' @param d1_pt First point of first line segment (Point object)
#' @param d2_pt Second point of first line segment (Point object)
#' @param m1_pos First point of second line segment (Point object)
#' @param m2_pos Second point of second line segment (Point object)
#' @return Point object representing intersection point or NULL if no intersection exists
find.line.intersection <- function(d1_pt, d2_pt, m1_pos, m2_pos) {
  # Calculate direction vectors
  dx1 <- d2_pt[1] - d1_pt[1]
  dy1 <- d2_pt[2] - d1_pt[2]
  dx2 <- m2_pos[1] - m1_pos[1]
  dy2 <- m2_pos[2] - m1_pos[2]

  # Calculate determinant to check if lines are parallel
  determinant <- dx1 * dy2 - dy1 * dx2

  # If determinant is close to zero, lines are parallel
  if (abs(determinant) < 1e-10) {
    return(NULL)
  }

  # Calculate differences between starting points
  dx3 <- d1_pt[1] - m1_pos[1]
  dy3 <- d1_pt[2] - m1_pos[2]

  # Calculate parameters for both lines
  t <- (dx2 * dy3 - dy2 * dx3) / determinant
  s <- (dx1 * dy3 - dy1 * dx3) / determinant

  # Check if intersection point lies within both segments
  if (!(t >= 0 && t <= 1 && s >= 0 && s <= 1)) {
    return(NULL)
  }

  # Calculate intersection point
  intersection_x <- d1_pt[1] + t * dx1
  intersection_y <- d1_pt[2] + t * dy1

  return(c(intersection_x, intersection_y))
}

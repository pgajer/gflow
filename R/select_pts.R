#' Select Points in 3D Space Using RGL
#'
#' A modified version of \code{rgl::selectpoints3d()} that provides more stable
#' point selection in 3D visualizations without crashing.
#'
#' @param objects A vector of rgl object ID values to search within.
#'   Defaults to all objects in the current rgl scene.
#' @param value Logical. If \code{TRUE} (default), returns the coordinates
#'   of selected points. If \code{FALSE}, returns their indices.
#' @param closest Logical. If \code{TRUE} (default), returns the points
#'   closest to the selection region when no points fall exactly within it.
#' @param multiple Logical or function. If \code{TRUE}, allows multiple
#'   selections. If a function, it should accept the current selection
#'   and return \code{TRUE} to continue selecting or \code{FALSE} to stop.
#' @param ... Additional parameters passed to \code{rgl::select3d()}.
#'
#' @return If \code{value = TRUE}, returns a matrix with columns 'x', 'y', 'z'
#'   containing the coordinates of selected points. If \code{value = FALSE},
#'   returns a matrix with columns 'id' and 'index' identifying the selected
#'   points within their respective objects.
#'
#' @details This function provides an interactive method for selecting points
#'   in 3D space. It improves upon the original \code{selectpoints3d()} by
#'   adding better error handling and preventing crashes that could occur
#'   with certain object configurations.
#'
#' @importFrom rgl ids3d rgl.attrib select3d rgl.user2window
#'
#' @examples
#' \dontrun{
#' library(rgl)
#' # Create a simple 3D scatter plot
#' x <- rnorm(100)
#' y <- rnorm(100)
#' z <- rnorm(100)
#' plot3d(x, y, z)
#'
#' # Select points interactively
#' selected <- points3d.select()
#'
#' # Highlight selected points
#' points3d(selected, col = "red", size = 10)
#' }
#'
#' @export
points3d.select <- function(objects = rgl::ids3d()$id,
                           value = TRUE,
                           closest = TRUE,
                           multiple = TRUE,
                           ...) {
    # Input validation
    if (!is.logical(value) || length(value) != 1) {
        stop("'value' must be a single logical value")
    }
    if (!is.logical(closest) || length(closest) != 1) {
        stop("'closest' must be a single logical value")
    }

    # Initialize result based on return type
    if (value) {
        result <- cbind(x = numeric(0), y = numeric(0), z = numeric(0))
    } else {
        result <- cbind(id = integer(0), index = integer(0))
    }

    rdist <- function(x) x  # Identity function for distance
    first <- TRUE

    # Main selection loop
    while (first || is.function(multiple) || isTRUE(multiple)) {
        # Get selection function
        f <- rgl::select3d(...)
        if (is.null(f)) {
            break
        }

        e <- environment(f)
        dist <- Inf
        prev <- nrow(result)

        # Process each object
        for (id in objects) {
            # Get vertices for current object
            verts <- try(rgl::rgl.attrib(id, "vertices"), silent = TRUE)
            if (inherits(verts, "try-error") || nrow(verts) < 1) {
                next
            }

            # Check for hits
            hits <- f(verts)
            if (any(hits)) {
                dist <- 0
            } else if (closest && dist > 0 && nrow(verts) > 0) {
                # Find closest points if no exact hits
                wincoords <- rgl::rgl.user2window(verts, projection = e$proj)
                wz <- wincoords[, 3]
                keep <- (0 <= wz) & (wz <= 1)
                wincoords <- wincoords[keep, , drop = FALSE]

                if (nrow(wincoords) == 0) {
                    next
                }

                wx <- wincoords[, 1]
                xdist <- ifelse(wx < e$llx, (wx - e$llx)^2,
                               ifelse(wx < e$urx, 0, (wx - e$urx)^2))

                wy <- wincoords[, 2]
                ydist <- ifelse(wy < e$lly, (wy - e$lly)^2,
                               ifelse(wy < e$ury, 0, (wy - e$ury)^2))

                dists <- xdist + ydist
                hits <- (dists < dist) & (dists == min(dists))
                dist <- min(c(dist, dists))
            }

            if (!any(hits)) {
                next
            }

            # Update results if closer points found
            if (prev > 0 && nrow(result) > prev && rdist(1) > dist) {
                result <- result[seq_len(prev), , drop = FALSE]
            }

            # Add selected points to results
            if (value) {
                result <- rbind(result, verts[hits, ])
            } else {
                result <- rbind(result, cbind(id = id, index = which(hits)))
            }

            # Check if multiple selection should continue
            if (is.function(multiple) && nrow(result) > prev) {
                new_selection <- result[(prev + 1):nrow(result), , drop = FALSE]
                if (!multiple(new_selection)) {
                    break
                }
            }

            rdist <- function(x) dist
            first <- FALSE
        }

        # Remove duplicates if returning coordinates
        if (value) {
            result <- unique(result)
        }
    }

    return(result)
}


#' Display Profile of Top Components
#'
#' Shows a profile of the columns with the largest mean values from a
#' selected row of a matrix, typically used for analyzing biomarker abundances.
#'
#' @param i Integer. The row index at which to extract the profile.
#' @param Z Matrix. A matrix of values (e.g., biomarker abundances) where
#'   rows represent samples and columns represent features.
#' @param n.comp Integer. The number of top components to include in the
#'   profile (default: 5).
#'
#' @return A named numeric vector containing the top \code{n.comp} features
#'   sorted by their mean values in descending order.
#'
#' @details This function extracts a single row from the input matrix and
#'   returns the features with the highest mean values. If the extracted
#'   row is a vector, it's converted to a single-row matrix for consistent
#'   processing.
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' Z <- matrix(rnorm(100, mean = 10, sd = 2), nrow = 10, ncol = 10)
#' colnames(Z) <- paste0("Feature", 1:10)
#'
#' # Get profile for row 1
#' profile <- show.profile(1, Z, n.comp = 3)
#' print(profile)
#'
#' @export
show.profile <- function(i, Z, n.comp = 5) {
    # Input validation
    if (!is.numeric(i) || length(i) != 1 || i < 1) {
        stop("'i' must be a single positive integer")
    }
    if (!is.matrix(Z) && !is.data.frame(Z)) {
        stop("'Z' must be a matrix or data frame")
    }
    if (i > nrow(Z)) {
        stop(sprintf("'i' (%d) exceeds number of rows in 'Z' (%d)", i, nrow(Z)))
    }
    if (!is.numeric(n.comp) || length(n.comp) != 1 || n.comp < 1) {
        stop("'n.comp' must be a single positive integer")
    }

    # Extract row
    Y <- Z[i, , drop = FALSE]

    # Calculate column means and sort
    col_means <- apply(Y, 2, mean, na.rm = TRUE)
    n.comp <- min(n.comp, length(col_means))

    # Get top components
    top_components <- sort(col_means, decreasing = TRUE)[seq_len(n.comp)]

    return(as.matrix(top_components, ncol = 1))
}


#' Interactive 3D Point Selection with Profile Display
#'
#' Allows interactive selection of points in a 3D visualization and optionally
#' displays profiles of the selected points based on associated data.
#'
#' @param X Numeric matrix with 3 columns representing x, y, z coordinates
#'   of points in 3D space.
#' @param Z Optional matrix or data frame where rows correspond to points in
#'   \code{X} and columns represent features (e.g., species abundances).
#' @param n.comp Integer. Number of top features to show in profiles (default: 5).
#' @param show.profiles Logical. Whether to print profiles to console (default: FALSE).
#' @param use.pts Logical. If TRUE, displays selected points as points;
#'   if FALSE, displays as spheres (default: FALSE).
#' @param color Character string specifying the color for highlighting
#'   selected points (default: "red").
#' @param alpha Numeric between 0 and 1. Transparency of spheres when
#'   \code{use.pts = FALSE} (default: 0.3).
#' @param size Numeric. Size of points when \code{use.pts = TRUE} (default: 0.1).
#' @param radius Numeric. Radius of spheres when \code{use.pts = FALSE} (default: 0.1).
#'
#' @return Invisibly returns a list with two components:
#'   \item{idx}{Integer vector of indices of selected points in \code{X}}
#'   \item{prof}{Matrix showing the mean profile of selected points' features,
#'     or NA if \code{Z} is not provided}
#'
#' @details This function provides an interactive interface for selecting
#'   points in an existing 3D rgl plot. Selected points are highlighted
#'   and their profiles can be computed from associated data.
#'
#' @importFrom rgl points3d spheres3d
#'
#' @examples
#' \dontrun{
#' library(rgl)
#' # Create example 3D data
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Z <- matrix(rpois(n * 10, lambda = 5), ncol = 10)
#' colnames(Z) <- paste0("Species", 1:10)
#'
#' # Create 3D plot
#' plot3d(X, col = "blue")
#'
#' # Interactively select points
#' selection <- pts3d.select(X, Z, show.profiles = TRUE)
#'
#' # Access results
#' selected_indices <- selection$idx
#' profile_data <- selection$prof
#' }
#'
#' @export
pts3d.select <- function(X,
                         Z = NULL,
                         n.comp = 5,
                         show.profiles = FALSE,
                         use.pts = FALSE,
                         color = "red",
                         alpha = 0.3,
                         size = 0.1,
                         radius = 0.1) {

    # Input validation
    if (!is.matrix(X) || ncol(X) != 3) {
        stop("'X' must be a matrix with exactly 3 columns")
    }
    if (!is.null(Z) && nrow(Z) != nrow(X)) {
        stop("'Z' must have the same number of rows as 'X'")
    }
    if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
        stop("'alpha' must be a numeric value between 0 and 1")
    }

    # Define selection callback function
    selection_callback <- function(ids) {
        indices <- ids[, "index"]
        points_to_highlight <- X[indices, , drop = FALSE]

        if (use.pts) {
            rgl::points3d(points_to_highlight, col = color, size = size)
        } else {
            rgl::spheres3d(points_to_highlight, col = color,
                          alpha = alpha, radius = radius)
        }
        return(TRUE)
    }

    # Perform interactive selection
    r <- points3d.select(value = FALSE, multiple = selection_callback)

    # Extract unique indices
    ii <- unique(r[, "index"])
    cat("Number of points selected:", length(ii), "\n")

    # Calculate profiles if data provided
    prof <- NA
    if (!is.null(Z) && length(ii) > 0) {
        Y <- Z[ii, , drop = FALSE]

        # Calculate mean profile
        col_means <- colMeans(Y, na.rm = TRUE)
        n.comp <- min(n.comp, length(col_means))
        prof <- sort(col_means, decreasing = TRUE)[seq_len(n.comp)]
        prof <- as.matrix(prof, ncol = 1)

        if (show.profiles) {
            cat("\nProfile of selected points:\n")
            print(prof)
        }
    }

    invisible(list(idx = ii, prof = prof))
}


#' Simple 3D Point Selection
#'
#' Provides a simplified interface for selecting points in 3D space without
#' visual feedback during selection.
#'
#' @param X Numeric matrix with 3 columns representing x, y, z coordinates
#'   of points in 3D space.
#'
#' @return Integer vector containing the indices of selected points.
#'
#' @details This function provides a minimal interface for point selection
#'   in 3D space. Unlike \code{pts3d.select}, it doesn't provide visual
#'   feedback during selection. It assumes the points in \code{X} correspond
#'   to vertices in the current rgl scene.
#'
#' @importFrom rgl ids3d rgl.attrib select3d
#'
#' @examples
#' \dontrun{
#' library(rgl)
#' # Create example 3D data
#' X <- matrix(rnorm(300), ncol = 3)
#'
#' # Plot points
#' plot3d(X, col = "blue")
#'
#' # Select points
#' selected_indices <- simple.pts3d.select(X)
#'
#' # Highlight selected points
#' points3d(X[selected_indices, ], col = "red", size = 10)
#' }
#'
#' @export
simple.pts3d.select <- function(X) {
    # Input validation
    if (!is.matrix(X) || ncol(X) != 3) {
        stop("'X' must be a matrix with exactly 3 columns")
    }

    # Get all object IDs in the scene
    ids <- rgl::ids3d()$id

    if (length(ids) == 0) {
        stop("No objects found in the current rgl scene")
    }

    # Find the object that matches our data
    matched_id <- NULL
    for (id in ids) {
        verts <- try(rgl::rgl.attrib(id, "vertices"), silent = TRUE)
        if (!inherits(verts, "try-error") &&
            nrow(verts) == nrow(X) &&
            ncol(verts) == 3) {
            matched_id <- id
            break
        }
    }

    if (is.null(matched_id)) {
        stop("Could not find an rgl object matching the dimensions of 'X'")
    }

    # Get vertices and perform selection
    verts <- rgl::rgl.attrib(matched_id, "vertices")
    f <- rgl::select3d(button = "left")

    if (is.null(f)) {
        warning("Selection was cancelled")
        return(integer(0))
    }

    hits <- f(verts)

    return(which(hits))
}

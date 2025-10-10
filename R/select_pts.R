#' Select 3D Points and Show Profiles
#'
#' Interactively select points in 3D space by drawing rectangles around them,
#' and optionally compute their profiles from associated data.
#'
#' @param X Numeric matrix with 3 columns representing 3D coordinates.
#'   **IMPORTANT**: Must be the exact same matrix used for plotting (e.g., in
#'   \code{plot3D.cont()} or \code{plot3D.plain()}) to ensure correct index mapping.
#' @param Z Optional data matrix with the same number of rows as \code{X} for computing profiles.
#' @param n.comp Number of top components to include in profiles (default: 5).
#' @param show.profiles Logical; whether to display profiles after selection completes.
#' @param use.pts Logical; if TRUE, use points for highlighting; otherwise use spheres.
#' @param color Color for highlighting selected points (default: "red").
#' @param alpha Transparency for highlighted spheres, 0-1 (default: 0.3).
#' @param size Size for highlighted points when \code{use.pts = TRUE} (default: 10).
#' @param radius Radius for highlighted spheres when \code{use.pts = FALSE} (default: 0.1).
#' @param allow_headless Logical; if TRUE, bypasses the headless/null-device check
#'   (selection still won't work without a real OpenGL device). Default FALSE.
#'
#' @return A list with:
#'   \item{idx}{Integer vector of selected row indices in \code{X}.}
#'   \item{prof}{1xk matrix of the mean profile over selected rows (top \code{n.comp} features),
#'               or \code{NA} if \code{Z} is not provided or nothing is selected.}
#'
#' @details
#' This function is **interactive** and requires a real rgl OpenGL window.
#' Selection is performed by **drawing rectangles** around points (click and drag),
#' not by single clicks. As you select points, they are highlighted with the specified
#' color and style. Right-click or press ESC when finished selecting.
#'
#' **Critical**: The \code{X} matrix must be identical to the matrix used in your
#' plotting command. For example, if you called \code{plot3D.cont(graph.3d, ...)},
#' you must use \code{X = graph.3d} here, not a different matrix.
#'
#' On CRAN/rhub or other headless environments (where \code{options(rgl.useNULL)=TRUE}),
#' interactive selection is unavailable and the function stops with an informative error.
#'
#' @examples
#' \dontrun{
#' if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
#'   set.seed(123)
#'   n <- 100
#'   X <- matrix(rnorm(n * 3), ncol = 3)
#'   Z <- matrix(rpois(n * 10, lambda = 5), ncol = 10)
#'   colnames(Z) <- paste0("Feature", 1:10)
#'
#'   # Plot the data
#'   rgl::open3d()
#'   rgl::plot3d(X, col = "blue")
#'
#'   # Select points and compute profiles - use the SAME X matrix!
#'   sel <- select3D.points.profiles(X, Z, show.profiles = TRUE)
#'
#'   # Access results
#'   print(sel$idx)  # Selected row indices
#'   print(sel$prof) # Mean profile of selected points
#' }
#' }
#' @export
select3D.points.profiles <- function(X,
                                     Z = NULL,
                                     n.comp = 5,
                                     show.profiles = FALSE,
                                     use.pts = FALSE,
                                     color = "red",
                                     alpha = 0.3,
                                     size = 10,
                                     radius = 0.1,
                                     allow_headless = FALSE) {

  # rgl gating
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl'. ",
         "Install with install.packages('rgl').", call. = FALSE)
  }

  # Input validation
  if (!is.matrix(X) || ncol(X) != 3L) {
    stop("'X' must be a numeric matrix with exactly 3 columns", call. = FALSE)
  }
  storage.mode(X) <- "double"
  if (!is.null(Z) && nrow(Z) != nrow(X)) {
    stop("'Z' must have the same number of rows as 'X'", call. = FALSE)
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a numeric value between 0 and 1", call. = FALSE)
  }

  use_null <- (!interactive()) ||
      identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
      (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
  old_opt <- options(rgl.useNULL = use_null)
  on.exit(options(old_opt), add = TRUE)

  # Require active device and scene contents
  if (rgl::rgl.cur() == 0) {
    stop("No active rgl device. Plot your 3D points first (e.g., rgl::plot3d(X)).",
         call. = FALSE)
  }
  ids_df <- rgl::ids3d()
  if (NROW(ids_df) == 0) {
    stop("No objects found in the current rgl scene.", call. = FALSE)
  }

  # Callback to highlight selected points as the user makes selections
  selection_callback <- function(ids) {
    if (is.null(ids) || NROW(ids) == 0) return(TRUE)
    sel_idx <- ids[, "index", drop = TRUE]
    if (!length(sel_idx)) return(TRUE)
    pts <- X[sel_idx, , drop = FALSE]
    if (isTRUE(use.pts)) {
      rgl::points3d(pts, col = color, size = size)
    } else {
      rgl::spheres3d(pts, col = color, alpha = alpha, radius = radius)
    }
    TRUE  # Continue selecting
  }

  # Perform interactive selection
  message("Draw rectangles around points to select them. Right-click when done.")
  r <- select3D.points(value = FALSE, multiple = selection_callback)

  if (is.null(r) || NROW(r) == 0) {
    message("No points selected.")
    return(list(idx = integer(0), prof = NA))
  }

  # Unique selected indices
  ii <- sort(unique(as.integer(r[, "index", drop = TRUE])))
  message(sprintf("Number of points selected: %d", length(ii)))

  # Compute profile if Z is provided and something was selected
  prof <- NA
  if (!is.null(Z) && length(ii) > 0) {
    Y <- Z[ii, , drop = FALSE]
    mean_prof <- colMeans(Y, na.rm = TRUE)
    ord <- order(abs(mean_prof), decreasing = TRUE)
    k <- min(n.comp, length(mean_prof))
    top <- ord[seq_len(k)]
    prof_vec <- mean_prof[top]
    # return a 1xk matrix as documented
    prof <- matrix(prof_vec, nrow = 1,
                   dimnames = list("mean", names(mean_prof)[top]))

    if (isTRUE(show.profiles)) {
      message("\nProfile of selected points (top ", k, " features):")
      print(round(prof, 3))
    }
  }

  invisible(list(idx = ii, prof = prof))
}

#' Select Points in 3D Space
#'
#' An interactive function for selecting points in 3D space by drawing rectangles
#' around them. This is a modified implementation that provides stable interaction
#' and proper handling of edge cases.
#'
#' @param objects Integer vector of rgl object IDs to consider for selection.
#'   Defaults to all objects in the current scene.
#' @param value Logical; if TRUE, return coordinates of selected points.
#'   If FALSE (default), return object IDs and indices.
#' @param closest Logical; currently not fully implemented. Reserved for future use.
#' @param multiple Logical or function; if TRUE, allow multiple selections.
#'   If a function, it's called after each selection with the newly added rows
#'   (as a matrix with columns matching the \code{value} parameter);
#'   return FALSE from the function to stop.
#' @param allow_headless Logical; if TRUE, bypasses the headless/null-device
#'   check (selection still won't work without a real OpenGL window). Default FALSE.
#'
#' @return If \code{value = TRUE}, a matrix with columns \code{x,y,z} of selected
#'   coordinates. If \code{value = FALSE}, an integer matrix with columns
#'   \code{id,index} identifying selected points, where \code{index} corresponds
#'   to the row number in the original plotted matrix.
#'
#' @details
#' This function is **interactive** and requires a real rgl OpenGL device.
#' Selection is performed by **drawing rectangles** around points (click and drag),
#' not by single clicks. Continue drawing rectangles to select multiple points,
#' and right-click or press ESC when finished.
#'
#' On CRAN/rhub or other headless environments (where \code{options(rgl.useNULL)=TRUE})
#' or in non-interactive sessions (\code{interactive()==FALSE}), the function stops
#' with an informative error.
#'
#' @examples
#' \dontrun{
#' if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(300), ncol = 3)
#'   rgl::open3d()
#'   rgl::plot3d(X)
#'
#'   # Return indices of selected points (default):
#'   result <- select3D.points()
#'   selected_rows <- result[, "index"]
#'
#'   # Or return coordinates of selected points:
#'   coords <- select3D.points(value = TRUE)
#' }
#' }
#' @export
select3D.points <- function(objects = NULL,
                            value = TRUE,
                            closest = TRUE,
                            multiple = TRUE,
                            allow_headless = FALSE) {

  # rgl gating
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl'. ",
         "Install with install.packages('rgl').", call. = FALSE)
  }

  use_null <- (!interactive()) ||
      identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
      (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
  old_opt <- options(rgl.useNULL = use_null)
  on.exit(options(old_opt), add = TRUE)

  # Check for active device
  if (rgl::rgl.cur() == 0) {
    stop("No active rgl device.", call. = FALSE)
  }

  # Default objects to all scene objects
  if (is.null(objects)) {
    ids_df <- rgl::ids3d()
    if (NROW(ids_df) == 0) stop("No objects found in the current rgl scene.", call. = FALSE)
    objects <- ids_df$id
  }

  message("Draw small rectangles around points to select them (click and drag).")
  message("Right-click or press ESC when done.")

  # Initialize result
  if (isTRUE(value)) {
    result <- matrix(nrow = 0, ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  } else {
    result <- matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c("id", "index")))
  }

  repeat {
    prev_n <- nrow(result)

    # Get selection rectangle
    f <- rgl::select3d()
    if (is.null(f)) break  # User cancelled

    for (id in objects) {
      verts <- try(rgl::rgl.attrib(id, "vertices"), silent = TRUE)
      if (inherits(verts, "try-error") || !is.matrix(verts) ||
          nrow(verts) == 0L || ncol(verts) < 3L) {
        next
      }

      # Test which vertices are inside selection
      hits <- try(f(verts), silent = TRUE)
      if (inherits(hits, "try-error") || length(hits) != nrow(verts)) {
        hits <- rep(FALSE, nrow(verts))
      }

      # CRITICAL: Handle NA values
      hits[is.na(hits)] <- FALSE

      if (any(hits)) {
        if (isTRUE(value)) {
          result <- rbind(result, verts[hits, 1:3, drop = FALSE])
        } else {
          result <- rbind(result, cbind(id = rep.int(id, sum(hits)),
                                       index = which(hits)))
        }
      }
    }

    # Call the callback function if provided
    if (is.function(multiple) && nrow(result) > prev_n) {
      new_rows <- result[(prev_n + 1):nrow(result), , drop = FALSE]
      cont <- try(multiple(new_rows), silent = TRUE)
      if (inherits(cont, "try-error") || isFALSE(cont)) {
        # Callback returned FALSE or error, stop selecting
        result <- result[seq_len(prev_n), , drop = FALSE]
        break
      }
    }

    if (!is.function(multiple) && !isTRUE(multiple)) break
  }

  # Ensure proper types
  if (!isTRUE(value) && nrow(result) > 0) {
    result[, "id"] <- as.integer(result[, "id"])
    result[, "index"] <- as.integer(result[, "index"])
  }

  message(sprintf("Selected %d point(s)", nrow(result)))
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

#' Select Points in 3D Space
#'
#' A modified version of \code{rgl::selectpoints3d()} that provides more stable
#' interaction for selecting points in 3D space.
#'
#' @param objects Integer vector of rgl object IDs to consider for selection.
#'   Defaults to all objects in the current scene.
#' @param value Logical; if TRUE, return coordinates of selected points.
#'   If FALSE, return object IDs and indices.
#' @param closest Logical; if TRUE, pick the closest point when clicking near
#'   but not directly on a point.
#' @param multiple Logical or function; if TRUE, allow multiple selections.
#'   If a function, it's called after each selection with the newly added rows;
#'   return FALSE from the function to stop.
#' @param allow_headless Logical; if TRUE, bypasses the headless/null-device
#'   check (selection still won't work without a real OpenGL window). Default FALSE.
#' @param ... Additional parameters passed to \code{rgl::select3d()}.
#'
#' @return If \code{value = TRUE}, a matrix with columns \code{x,y,z} of selected
#'   coordinates. If \code{value = FALSE}, an integer matrix with columns
#'   \code{id,index} identifying selected points.
#'
#' @details
#' This function is **interactive** and requires a real rgl OpenGL device.
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
#'   # Return coordinates of selected points:
#'   coords <- points3D.select()
#'   # Or return (id, index) pairs:
#'   # idxmat <- select3D.points(value = FALSE)
#' }
#' }
#' @export
select3D.points <- function(objects = NULL,
                            value = TRUE,
                            closest = TRUE,
                            multiple = TRUE,
                            allow_headless = FALSE,
                            ...) {

  # rgl gating
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl'. ",
         "Install with install.packages('rgl').", call. = FALSE)
  }
  # Interactive selection cannot run headless or in non-interactive sessions
  if (!isTRUE(allow_headless) && isTRUE(getOption("rgl.useNULL", FALSE))) {
    stop("Interactive selection is unavailable with the rgl null (off-screen) device. ",
         "Run in an interactive session with a real rgl window ",
         "(options(rgl.useNULL)=FALSE).", call. = FALSE)
  }
  if (!interactive()) {
    stop("Interactive selection requires an interactive R session (interactive() == TRUE).",
         call. = FALSE)
  }
  if (rgl::rgl.cur() == 0) {
    stop("No active rgl device. Plot your 3D points first (e.g., rgl::plot3d(...)).",
         call. = FALSE)
  }

  # Default objects to all scene objects
  if (is.null(objects)) {
    ids_df <- rgl::ids3d()
    if (NROW(ids_df) == 0) stop("No objects found in the current rgl scene.", call. = FALSE)
    objects <- ids_df$id
  }
  # Coerce/validate object IDs
  objects <- as.integer(objects)
  objects <- objects[is.finite(objects)]
  if (!length(objects)) stop("'objects' must contain valid rgl object IDs.", call. = FALSE)

  # Validate flags
  if (!is.logical(value)   || length(value)   != 1L) stop("'value' must be a single logical")
  if (!is.logical(closest) || length(closest) != 1L) stop("'closest' must be a single logical")

  # Initialize result
  if (isTRUE(value)) {
    result <- cbind(x = numeric(0), y = numeric(0), z = numeric(0))
  } else {
    result <- cbind(id = integer(0), index = integer(0))
  }

  first <- TRUE
  repeat {
    # Obtain selection function (user draws rectangle/region)
    f <- rgl::select3d(...)
    if (is.null(f)) break  # user canceled/finished

    env <- environment(f)
    have_proj <- is.environment(env) && all(c("proj", "llx", "lly") %in% ls(env, all.names = TRUE))

    prev_n <- nrow(result)
    picked_any <- FALSE
    best_dist <- Inf

    for (id in objects) {
      verts <- try(rgl::rgl.attrib(id, "vertices"), silent = TRUE)
      if (inherits(verts, "try-error") || !is.matrix(verts) || nrow(verts) == 0L || ncol(verts) < 3L)
        next

      # Direct hits via selection function
      hits <- try(f(verts), silent = TRUE)
      if (inherits(hits, "try-error") || length(hits) != nrow(verts)) hits <- rep(FALSE, nrow(verts))

      if (any(hits)) {
        picked_any <- TRUE
        if (isTRUE(value)) {
          result <- rbind(result, verts[hits, 1:3, drop = FALSE])
        } else {
          result <- rbind(result, cbind(id = rep.int(id, sum(hits)), index = which(hits)))
        }
        next
      }

      # If no hits, optionally choose the closest point to the selection
      if (isTRUE(closest) && have_proj && is.finite(best_dist)) {
        wc <- try(rgl::rgl.user2window(verts, projection = env$proj), silent = TRUE)
        if (!inherits(wc, "try-error") && is.matrix(wc) && ncol(wc) >= 3L) {
          keep <- (wc[, 3] >= 0) & (wc[, 3] <= 1)
          if (any(keep)) {
            wx <- wc[keep, 1] - env$llx
            wy <- wc[keep, 2] - env$lly
            d  <- sqrt(wx * wx + wy * wy)
            if (length(d) && min(d) < best_dist) {
              best_dist <- min(d)
              imin <- which.min(d)
              if (isTRUE(value)) {
                result <- rbind(result, verts[keep, 1:3, drop = FALSE][imin, , drop = FALSE])
              } else {
                idx_keep <- which(keep)
                result <- rbind(result, cbind(id = id, index = idx_keep[imin]))
              }
              picked_any <- TRUE
            }
          }
        }
      }
    } # for each object

    # Invoke callback for multiple=function
    if (is.function(multiple) && nrow(result) > prev_n) {
      new_rows <- result[(prev_n + 1):nrow(result), , drop = FALSE]
      cont <- try(multiple(new_rows), silent = TRUE)
      if (inherits(cont, "try-error") || isFALSE(cont)) {
        result <- result[seq_len(prev_n), , drop = FALSE]
        break
      }
    }

    # Single selection mode: stop after first pick
    if (!is.function(multiple) && !isTRUE(multiple)) break

    # If user didn’t pick anything this round, allow trying again; otherwise loop
    first <- FALSE
  }

  # Ensure integer columns for id/index form
  if (!isTRUE(value) && nrow(result)) {
    result[, "id"]    <- as.integer(result[, "id"])
    result[, "index"] <- as.integer(result[, "index"])
  }

  result
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

#' Select 3D Points and Show Profiles
#'
#' Interactively select points in 3D space and optionally compute their profiles
#' from associated data.
#'
#' @param X Numeric matrix with 3 columns representing 3D coordinates.
#' @param Z Optional data matrix with the same number of rows as \code{X} for computing profiles.
#' @param n.comp Number of top components to include in profiles (default: 5).
#' @param show.profiles Logical; whether to display profiles during selection.
#' @param use.pts Logical; if TRUE, use points for highlighting; otherwise use spheres.
#' @param color Color for highlighting selected points.
#' @param alpha Transparency for highlighted spheres (0–1).
#' @param size Size for highlighted points (when \code{use.pts = TRUE}).
#' @param radius Radius for highlighted spheres (when \code{use.pts = FALSE}).
#' @param allow_headless Logical; if TRUE, bypasses the headless/null-device check
#'   (selection still won't work without a real OpenGL device). Default FALSE.
#'
#' @return A list with:
#'   \item{idx}{Integer vector of selected row indices in \code{X}.}
#'   \item{prof}{1×k matrix of the mean profile over selected rows (top \code{n.comp} features),
#'               or \code{NA} if \code{Z} is not provided or nothing is selected.}
#'
#' @details
#' This function is **interactive** and requires a real rgl OpenGL window. On CRAN/rhub
#' or other headless environments (where \code{options(rgl.useNULL)=TRUE}), interactive
#' selection is unavailable and the function stops with an informative error.
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
#'   rgl::open3d()
#'   rgl::plot3d(X, col = "blue")
#'
#'   sel <- select3D.points.profiles(X, Z, show.profiles = TRUE)
#'   rgl::points3d(X[sel$idx, , drop = FALSE], col = "red", size = 10)
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
                                     size = 0.1,
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

  # Interactive selection cannot work with the null (off-screen) device
  if (!isTRUE(allow_headless) && isTRUE(getOption("rgl.useNULL", FALSE))) {
    stop("Interactive selection is unavailable with the rgl null (off-screen) device. ",
         "Run in an interactive session with a real rgl window ",
         "(options(rgl.useNULL)=FALSE).", call. = FALSE)
  }

  # Require active device and scene contents
  if (rgl::rgl.cur() == 0) {
    stop("No active rgl device. Plot your 3D points first (e.g., rgl::plot3d(X)).",
         call. = FALSE)
  }
  ids_df <- rgl::ids3d()
  if (NROW(ids_df) == 0) {
    stop("No objects found in the current rgl scene.", call. = FALSE)
  }

  # Ensure selection helper exists
  if (!exists("points3D.select", mode = "function")) {
    stop("Required helper 'points3D.select()' not found in the package namespace.",
         call. = FALSE)
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
    TRUE
  }

  # Perform interactive selection (blocks until user finishes; API provided by your helper)
  r <- points3D.select(value = FALSE, multiple = selection_callback)

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
    mean_prof <- colMeans(Y)
    ord <- order(abs(mean_prof), decreasing = TRUE)
    k <- min(n.comp, length(mean_prof))
    top <- ord[seq_len(k)]
    prof_vec <- mean_prof[top]
    # return a 1×k matrix as documented
    prof <- matrix(prof_vec, nrow = 1,
                   dimnames = list("mean", colnames(Z)[top]))

    if (isTRUE(show.profiles)) {
      message("\nProfile of selected points (top ", k, " features):")
      print(round(prof, 3))
      if (exists("show.profile", mode = "function")) {
        for (idx in ii) {
          try(show.profile(idx, Z, n.comp), silent = TRUE)
        }
      }
    }
  }

  invisible(list(idx = ii, prof = prof))
}

#' Select Points in 3D Space
#'
#' A modified version of \code{rgl::selectpoints3d()} that provides more stable
#' interaction for selecting points in 3D space.
#'
#' @param objects Integer vector of rgl object IDs to consider for selection.
#'   Defaults to all objects in the current scene.
#' @param value Logical; if TRUE, return coordinates of selected points.
#'   If FALSE, return object IDs and indices.
#' @param closest Logical; if TRUE, pick the closest point when clicking near
#'   but not directly on a point.
#' @param multiple Logical or function; if TRUE, allow multiple selections.
#'   If a function, it's called after each selection with the newly added rows;
#'   return FALSE from the function to stop.
#' @param allow_headless Logical; if TRUE, bypasses the headless/null-device
#'   check (selection still won't work without a real OpenGL window). Default FALSE.
#' @param ... Additional parameters passed to \code{rgl::select3d()}.
#'
#' @return If \code{value = TRUE}, a matrix with columns \code{x,y,z} of selected
#'   coordinates. If \code{value = FALSE}, an integer matrix with columns
#'   \code{id,index} identifying selected points.
#'
#' @details
#' This function is **interactive** and requires a real rgl OpenGL device.
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
#'   # Return coordinates of selected points:
#'   coords <- points3D.select()
#'   # Or return (id, index) pairs:
#'   # idxmat <- points3d.select(value = FALSE)
#' }
#' }
#' @export
points3D.select <- function(objects = NULL,
                            value = TRUE,
                            closest = TRUE,
                            multiple = TRUE,
                            allow_headless = FALSE,
                            ...) {

  # rgl gating
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl'. ",
         "Install with install.packages('rgl').", call. = FALSE)
  }
  # Interactive selection cannot run headless or in non-interactive sessions
  if (!isTRUE(allow_headless) && isTRUE(getOption("rgl.useNULL", FALSE))) {
    stop("Interactive selection is unavailable with the rgl null (off-screen) device. ",
         "Run in an interactive session with a real rgl window ",
         "(options(rgl.useNULL)=FALSE).", call. = FALSE)
  }
  if (!interactive()) {
    stop("Interactive selection requires an interactive R session (interactive() == TRUE).",
         call. = FALSE)
  }
  if (rgl::rgl.cur() == 0) {
    stop("No active rgl device. Plot your 3D points first (e.g., rgl::plot3d(...)).",
         call. = FALSE)
  }

  # Default objects to all scene objects
  if (is.null(objects)) {
    ids_df <- rgl::ids3d()
    if (NROW(ids_df) == 0) stop("No objects found in the current rgl scene.", call. = FALSE)
    objects <- ids_df$id
  }
  # Coerce/validate object IDs
  objects <- as.integer(objects)
  objects <- objects[is.finite(objects)]
  if (!length(objects)) stop("'objects' must contain valid rgl object IDs.", call. = FALSE)

  # Validate flags
  if (!is.logical(value)   || length(value)   != 1L) stop("'value' must be a single logical")
  if (!is.logical(closest) || length(closest) != 1L) stop("'closest' must be a single logical")

  # Initialize result
  if (isTRUE(value)) {
    result <- cbind(x = numeric(0), y = numeric(0), z = numeric(0))
  } else {
    result <- cbind(id = integer(0), index = integer(0))
  }

  first <- TRUE
  repeat {
    # Obtain selection function (user draws rectangle/region)
    f <- rgl::select3d(...)
    if (is.null(f)) break  # user canceled/finished

    env <- environment(f)
    have_proj <- is.environment(env) && all(c("proj", "llx", "lly") %in% ls(env, all.names = TRUE))

    prev_n <- nrow(result)
    picked_any <- FALSE
    best_dist <- Inf

    for (id in objects) {
      verts <- try(rgl::rgl.attrib(id, "vertices"), silent = TRUE)
      if (inherits(verts, "try-error") || !is.matrix(verts) || nrow(verts) == 0L || ncol(verts) < 3L)
        next

      # Direct hits via selection function
      hits <- try(f(verts), silent = TRUE)
      if (inherits(hits, "try-error") || length(hits) != nrow(verts)) hits <- rep(FALSE, nrow(verts))

      if (any(hits)) {
        picked_any <- TRUE
        if (isTRUE(value)) {
          result <- rbind(result, verts[hits, 1:3, drop = FALSE])
        } else {
          result <- rbind(result, cbind(id = rep.int(id, sum(hits)), index = which(hits)))
        }
        next
      }

      # If no hits, optionally choose the closest point to the selection
      if (isTRUE(closest) && have_proj && is.finite(best_dist)) {
        wc <- try(rgl::rgl.user2window(verts, projection = env$proj), silent = TRUE)
        if (!inherits(wc, "try-error") && is.matrix(wc) && ncol(wc) >= 3L) {
          keep <- (wc[, 3] >= 0) & (wc[, 3] <= 1)
          if (any(keep)) {
            wx <- wc[keep, 1] - env$llx
            wy <- wc[keep, 2] - env$lly
            d  <- sqrt(wx * wx + wy * wy)
            if (length(d) && min(d) < best_dist) {
              best_dist <- min(d)
              imin <- which.min(d)
              if (isTRUE(value)) {
                result <- rbind(result, verts[keep, 1:3, drop = FALSE][imin, , drop = FALSE])
              } else {
                idx_keep <- which(keep)
                result <- rbind(result, cbind(id = id, index = idx_keep[imin]))
              }
              picked_any <- TRUE
            }
          }
        }
      }
    } # for each object

    # Invoke callback for multiple=function
    if (is.function(multiple) && nrow(result) > prev_n) {
      new_rows <- result[(prev_n + 1):nrow(result), , drop = FALSE]
      cont <- try(multiple(new_rows), silent = TRUE)
      if (inherits(cont, "try-error") || isFALSE(cont)) {
        result <- result[seq_len(prev_n), , drop = FALSE]
        break
      }
    }

    # Single selection mode: stop after first pick
    if (!is.function(multiple) && !isTRUE(multiple)) break

    # If user didn’t pick anything this round, allow trying again; otherwise loop
    first <- FALSE
  }

  # Ensure integer columns for id/index form
  if (!isTRUE(value) && nrow(result)) {
    result[, "id"]    <- as.integer(result[, "id"])
    result[, "index"] <- as.integer(result[, "index"])
  }

  result
}

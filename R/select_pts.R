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

#' Select Points in 3D from an HTML/WebGL View
#'
#' Interactively select points by drawing rectangular brushes in a browser-based
#' 3D view powered by \code{rglwidget}.
#'
#' @param X Numeric matrix/data.frame with exactly 3 columns representing 3D coordinates.
#' @param plot.type Plot style used to render the HTML widget:
#'   \code{"plain"}, \code{"cont"}, or \code{"cltrs"}.
#' @param y Optional numeric vector used when \code{plot.type = "cont"}.
#' @param cltr Optional cluster labels used when \code{plot.type = "cltrs"}.
#' @param plot.args Named list of additional arguments passed to the selected
#'   widget plotting helper (\code{plot3D.plain.widget()},
#'   \code{plot3D.cont.widget()}, or \code{plot3D.cltrs.widget()}).
#' @param title Title shown in the selection app.
#' @param widget.width,widget.height Width/height (pixels) of the 3D widget.
#' @param launch.browser Logical; passed to \code{shiny::runApp()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{idx}: integer indices of selected rows in \code{X}.
#'   \item \code{coords}: selected coordinates (\code{X[idx, ]}).
#'   \item \code{n.selected}: number of selected rows.
#'   \item \code{plot.type}: effective plot type used in the app.
#'   \item \code{cancelled}: logical flag indicating whether the app was cancelled.
#' }
#'
#' @details
#' This function is intended for environments where native \code{rgl} OpenGL
#' selection is unavailable, but browser/WebGL rendering works.
#'
#' Workflow:
#' \enumerate{
#'   \item Switch mouse mode to \code{"selecting"}.
#'   \item Draw a rectangle on the 3D view.
#'   \item Click \code{"Add Brushed Points"} to accumulate selection.
#'   \item Repeat as needed, then click \code{"Done"}.
#' }
#'
#' @examples
#' \dontrun{
#' if (interactive() &&
#'     requireNamespace("rgl", quietly = TRUE) &&
#'     requireNamespace("shiny", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(600), ncol = 3)
#'   out <- select3D.points.html(X, plot.type = "plain")
#'   print(out$idx)
#' }
#' }
#' @export
select3D.points.html <- function(X,
                                 plot.type = c("plain", "cont", "cltrs"),
                                 y = NULL,
                                 cltr = NULL,
                                 plot.args = list(),
                                 title = "Select 3D Points (HTML)",
                                 widget.width = 1700L,
                                 widget.height = 1000L,
                                 launch.browser = TRUE) {
  if (!interactive() && !isTRUE(launch.browser)) {
    stop(
      "Non-interactive session detected. Set launch.browser = TRUE ",
      "or run from an interactive R session.",
      call. = FALSE
    )
  }
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl'. ",
         "Install with install.packages('rgl').", call. = FALSE)
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("This function requires the optional package 'shiny'. ",
         "Install with install.packages('shiny').", call. = FALSE)
  }

  if (!is.matrix(X) && !is.data.frame(X)) stop("'X' must be a matrix or data frame.", call. = FALSE)
  X <- as.matrix(X)
  if (ncol(X) != 3L) stop("'X' must have exactly 3 columns.", call. = FALSE)
  storage.mode(X) <- "double"
  if (is.null(colnames(X))) colnames(X) <- c("x", "y", "z")

  plot.type <- match.arg(plot.type)
  if (!is.list(plot.args)) stop("'plot.args' must be a named list.", call. = FALSE)

  if (identical(plot.type, "cont")) {
    if (is.null(y)) stop("'y' must be provided when plot.type = 'cont'.", call. = FALSE)
    y <- as.numeric(y)
    if (length(y) != nrow(X)) stop("length('y') must match nrow(X).", call. = FALSE)
  }
  if (identical(plot.type, "cltrs") && !is.null(cltr) && length(cltr) != nrow(X)) {
    stop("length('cltr') must match nrow(X).", call. = FALSE)
  }

  widget.width <- as.integer(widget.width)
  widget.height <- as.integer(widget.height)
  if (!is.finite(widget.width) || widget.width <= 0L) stop("'widget.width' must be positive.", call. = FALSE)
  if (!is.finite(widget.height) || widget.height <= 0L) stop("'widget.height' must be positive.", call. = FALSE)

  brush.id <- "gflow_select3d_brush"
  scene.id <- "gflow_select3d_widget"
  force.selecting.layer <- function(ctx = NULL) {
    # Ensure left-button rectangular selection is active in the widget.
    rgl::par3d(mouseMode = c("none", "selecting", "zoom", "fov", "pull"))
    invisible(NULL)
  }

  build.widget <- function() {
    args <- plot.args

    # Enforce app-local rendering and dimensions.
    args$output.file <- NULL
    args$open.browser <- FALSE
    args$widget.width <- widget.width
    args$widget.height <- widget.height
    args$shiny.brush <- brush.id
    # In Shiny render calls prependContent() is ignored; disable HTML legends
    # by default to avoid noisy warnings.
    if (identical(plot.type, "cont") && is.null(args$legend.show)) {
      args$legend.show <- FALSE
    }
    if (identical(plot.type, "cltrs") && is.null(args$show.legend)) {
      args$show.legend <- FALSE
    }
    if (is.null(args$post.layers)) {
      args$post.layers <- list(force.selecting.layer)
    } else if (is.function(args$post.layers)) {
      args$post.layers <- list(force.selecting.layer, args$post.layers)
    } else if (is.list(args$post.layers)) {
      if (!is.null(args$post.layers$fun)) {
        args$post.layers <- list(force.selecting.layer, args$post.layers)
      } else {
        args$post.layers <- c(list(force.selecting.layer), args$post.layers)
      }
    } else {
      stop("'plot.args$post.layers' must be NULL, function, or list.", call. = FALSE)
    }

    w <- if (identical(plot.type, "plain")) {
      do.call(plot3D.plain.widget, c(list(X = X), args))
    } else if (identical(plot.type, "cont")) {
      do.call(plot3D.cont.widget, c(list(X = X, y = y), args))
    } else {
      do.call(plot3D.cltrs.widget, c(list(X = X, cltr = cltr), args))
    }
    # Keep a stable id so Shiny selection plumbing references the right widget.
    w$elementId <- scene.id
    # Defensive patch for some rgl/Shiny/browser combinations where
    # mouse handlers briefly see undefined subscene ids.
    if (requireNamespace("htmlwidgets", quietly = TRUE)) {
      w <- htmlwidgets::onRender(
        w,
        "function(el, x) {
           if (!this || typeof this.getObj !== 'function') return;
           if (this.__gflow_sel_patch__) return;
           var self = this;
           var oldGetObj = self.getObj;
           self.getObj = function(id) {
             if (typeof id === 'undefined' || id === null) {
               if (self.scene && typeof self.scene.rootSubscene !== 'undefined') {
                 id = self.scene.rootSubscene;
               }
             }
             return oldGetObj.call(self, id);
           };
           var oldWhich = self.whichSubscene;
           if (typeof oldWhich === 'function') {
             self.whichSubscene = function(coords) {
               var id = oldWhich.call(self, coords);
               if (typeof id === 'undefined' || id === null) {
                 if (self.scene && typeof self.scene.rootSubscene !== 'undefined') {
                   return self.scene.rootSubscene;
                 }
               }
               return id;
             };
           }
           self.__gflow_sel_patch__ = true;
         }"
      )
    }
    w
  }

  brush.to.indices <- function(brush) {
    if (is.null(brush) || length(brush) == 0L) return(integer(0))
    if (!is.null(brush$state) && identical(brush$state, "inactive")) return(integer(0))

    sel.fn <- try(rgl::selectionFunction3d(brush), silent = TRUE)
    if (inherits(sel.fn, "try-error")) return(integer(0))

    keep <- try(sel.fn(X), silent = TRUE)
    if (inherits(keep, "try-error") || length(keep) != nrow(X)) return(integer(0))
    keep <- as.logical(keep)
    keep[is.na(keep)] <- FALSE
    which(keep)
  }

  empty.result <- list(
    idx = integer(0),
    coords = X[integer(0), , drop = FALSE],
    n.selected = 0L,
    plot.type = plot.type,
    cancelled = TRUE
  )

  app <- shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::tags$h3(title),
      shiny::fluidRow(
        shiny::column(
          width = 9,
          shiny::tags$div(
            style = "margin-bottom:8px;color:#444;",
            "Selection mode is fixed: left-drag on the plot to draw a selection rectangle."
          ),
          rgl::rglwidgetOutput(
            scene.id,
            width = paste0(widget.width, "px"),
            height = paste0(widget.height, "px")
          )
        ),
        shiny::column(
          width = 3,
          shiny::actionButton("add", "Add Brushed Points"),
          shiny::tags$br(),
          shiny::tags$br(),
          shiny::actionButton("clear", "Clear Selection"),
          shiny::tags$br(),
          shiny::tags$br(),
          shiny::actionButton("done", "Done"),
          shiny::tags$span(" "),
          shiny::actionButton("cancel", "Cancel"),
          shiny::tags$hr(),
          shiny::verbatimTextOutput("status"),
          shiny::tableOutput("preview")
        )
      )
    ),
    server = function(input, output, session) {
      selected.idx <- shiny::reactiveVal(integer(0))

      reset.brush <- function() {
        if ("shinyResetBrush" %in% getNamespaceExports("rgl")) {
          try(rgl::shinyResetBrush(session, brush.id), silent = TRUE)
        }
      }

      brushed.idx <- shiny::reactive({
        brush.to.indices(input[[brush.id]])
      })

      output[[scene.id]] <- rgl::renderRglwidget({
        build.widget()
      })

      shiny::observeEvent(input$add, {
        idx <- brushed.idx()
        if (length(idx) > 0L) {
          merged <- sort(unique(c(selected.idx(), as.integer(idx))))
          selected.idx(merged)
        }
        reset.brush()
      })

      shiny::observeEvent(input$clear, {
        selected.idx(integer(0))
        reset.brush()
      })

      output$status <- shiny::renderText({
        cur <- selected.idx()
        now <- brushed.idx()
        paste0(
          "Mode hint: left-drag to draw selection rectangle\n",
          "Brushed now: ", length(now), " point(s)\n",
          "Selected total: ", length(cur), " point(s)\n",
          if (length(cur) > 0L) {
            paste0("Index range: [", min(cur), ", ", max(cur), "]")
          } else {
            "Index range: NA"
          }
        )
      })

      output$preview <- shiny::renderTable({
        idx <- selected.idx()
        if (length(idx) < 1L) return(NULL)
        show.idx <- head(idx, 20L)
        out <- data.frame(
          index = as.integer(show.idx),
          rowname = if (!is.null(rownames(X))) as.character(rownames(X)[show.idx]) else NA_character_,
          stringsAsFactors = FALSE
        )
        out
      }, striped = TRUE, hover = TRUE)

      shiny::observeEvent(input$done, {
        idx <- sort(unique(as.integer(selected.idx())))
        shiny::stopApp(list(
          idx = idx,
          coords = X[idx, , drop = FALSE],
          n.selected = length(idx),
          plot.type = plot.type,
          cancelled = FALSE
        ))
      })

      shiny::observeEvent(input$cancel, {
        shiny::stopApp(empty.result)
      })
    }
  )

  result <- shiny::runApp(app, launch.browser = isTRUE(launch.browser))
  if (is.null(result)) result <- empty.result
  if (is.null(result$idx)) result$idx <- integer(0)
  result$idx <- sort(unique(as.integer(result$idx)))
  if (is.null(result$coords)) {
    result$coords <- X[result$idx, , drop = FALSE]
  }
  if (is.null(result$n.selected)) result$n.selected <- as.integer(length(result$idx))
  if (is.null(result$plot.type)) result$plot.type <- plot.type
  if (is.null(result$cancelled)) result$cancelled <- FALSE

  invisible(result)
}

#' Select 3D Points Using a Plotly Browser View
#'
#' Interactively select points in a browser by clicking points in a 3D Plotly
#' scatter plot. This avoids the \code{rgl} JavaScript mouse-selection path.
#'
#' @param X Numeric matrix/data.frame with exactly 3 columns representing 3D coordinates.
#' @param plot.type Plot style used for color mapping:
#'   \code{"plain"}, \code{"cont"}, or \code{"cltrs"}.
#' @param y Optional numeric vector used when \code{plot.type = "cont"}.
#' @param cltr Optional cluster labels used when \code{plot.type = "cltrs"}.
#' @param subset Optional logical vector; points with \code{FALSE} are drawn in gray
#'   and lower opacity.
#' @param point.size Marker size.
#' @param title Title shown in the selection app.
#' @param end.labels Optional named character vector for endpoint labels.
#'   Names are row indices (1-based, relative to \code{X}).
#' @param launch.browser Logical; passed to \code{shiny::runApp()}.
#' @param save.dir Default directory used when saving selections.
#' @param save.prefix Filename prefix used for saved selection files.
#' @param save.on.done Logical; if \code{TRUE}, \code{Done} saves selection
#'   before closing the app.
#' @param allow.edit.save.dir Logical; if \code{TRUE}, the app shows an editable
#'   save-directory input box.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{idx}: integer indices of selected rows in \code{X}.
#'   \item \code{coords}: selected coordinates (\code{X[idx, ]}).
#'   \item \code{n.selected}: number of selected rows.
#'   \item \code{plot.type}: effective plot type.
#'   \item \code{cancelled}: logical flag indicating whether app was cancelled.
#'   \item \code{saved.files}: list with paths of last saved CSV/RDS (or \code{NULL}).
#' }
#'
#' @details
#' In the app, click a point to toggle selection. Selected points are highlighted
#' in red. Use \code{Undo}, \code{Clear}, \code{Save Now}, and \code{Done} controls.
#'
#' @examples
#' \dontrun{
#' if (interactive() &&
#'     requireNamespace("shiny", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'   X <- matrix(rnorm(600), ncol = 3)
#'   out <- select3D.points.plotly(X, plot.type = "plain")
#'   print(out$idx)
#' }
#' }
#' @export
select3D.points.plotly <- function(X,
                                   plot.type = c("plain", "cont", "cltrs"),
                                   y = NULL,
                                   cltr = NULL,
                                   subset = NULL,
                                   point.size = 3,
                                   title = "Select 3D Points (Plotly)",
                                   end.labels = NULL,
                                   launch.browser = TRUE,
                                   save.dir = tempdir(),
                                   save.prefix = "select3d_points",
                                   save.on.done = TRUE,
                                   allow.edit.save.dir = TRUE) {
  if (!interactive() && !isTRUE(launch.browser)) {
    stop(
      "Non-interactive session detected. Set launch.browser = TRUE ",
      "or run from an interactive R session.",
      call. = FALSE
    )
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("This function requires the optional package 'shiny'. ",
         "Install with install.packages('shiny').", call. = FALSE)
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("This function requires the optional package 'plotly'. ",
         "Install with install.packages('plotly').", call. = FALSE)
  }

  if (!is.matrix(X) && !is.data.frame(X)) stop("'X' must be a matrix or data frame.", call. = FALSE)
  X <- as.matrix(X)
  if (ncol(X) != 3L) stop("'X' must have exactly 3 columns.", call. = FALSE)
  storage.mode(X) <- "double"
  if (is.null(colnames(X))) colnames(X) <- c("x", "y", "z")

  plot.type <- match.arg(plot.type)
  n <- nrow(X)

  if (is.null(subset)) {
    subset <- rep(TRUE, n)
  } else {
    subset <- as.logical(subset)
    if (length(subset) != n) stop("length('subset') must match nrow(X).", call. = FALSE)
    subset[is.na(subset)] <- FALSE
  }

  if (!is.numeric(point.size) || length(point.size) != 1L || !is.finite(point.size) || point.size <= 0) {
    stop("'point.size' must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.character(save.dir) || length(save.dir) != 1L || !nzchar(save.dir)) {
    stop("'save.dir' must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.character(save.prefix) || length(save.prefix) != 1L || !nzchar(save.prefix)) {
    stop("'save.prefix' must be a non-empty character scalar.", call. = FALSE)
  }

  colors <- rep("gray45", n)
  legend.info <- NULL

  if (identical(plot.type, "cont")) {
    if (is.null(y)) stop("'y' must be provided when plot.type = 'cont'.", call. = FALSE)
    y <- as.numeric(y)
    if (length(y) != n) stop("length('y') must match nrow(X).", call. = FALSE)
    q <- quantize.cont.var(
      y,
      method = "quantile",
      n.levels = 11,
      na.color = "gray80"
    )
    colors <- unname(q$x.col.tbl[as.character(q$x.cat)])
    colors[is.na(colors)] <- "gray80"
    legend.info <- q$x.col.tbl
  } else if (identical(plot.type, "cltrs")) {
    if (is.null(cltr)) stop("'cltr' must be provided when plot.type = 'cltrs'.", call. = FALSE)
    cltr <- as.character(cltr)
    if (length(cltr) != n) stop("length('cltr') must match nrow(X).", call. = FALSE)
    cltr.labs <- unique(cltr)
    cltr.labs <- cltr.labs[order(cltr.labs)]
    pal <- grDevices::hcl.colors(length(cltr.labs), palette = "Dark 3")
    tbl <- setNames(pal, cltr.labs)
    colors <- unname(tbl[cltr])
    colors[is.na(colors)] <- "gray80"
    legend.info <- tbl
  }

  all.idx <- seq_len(n)
  idx.main <- all.idx[subset]
  idx.gray <- all.idx[!subset]

  if (is.null(end.labels)) {
    end.labels <- character(0)
  } else {
    end.labels <- as.character(end.labels)
    names(end.labels) <- names(end.labels)
    nm <- suppressWarnings(as.integer(names(end.labels)))
    ok <- is.finite(nm) & nm >= 1L & nm <= n
    end.labels <- end.labels[ok]
    names(end.labels) <- as.character(nm[ok])
  }

  empty.result <- list(
    idx = integer(0),
    coords = X[integer(0), , drop = FALSE],
    n.selected = 0L,
    plot.type = plot.type,
    cancelled = TRUE,
    saved.files = NULL
  )

  source.id <- "gflow_select3d_plotly_source"

  stamp.now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
  sanitize.prefix <- function(x) {
    y <- gsub("[^0-9A-Za-z_\\-]", "_", as.character(x))
    y <- gsub("_+", "_", y)
    y <- gsub("^_+|_+$", "", y)
    if (!nzchar(y)) y <- "select3d_points"
    y
  }
  write.selection.files <- function(idx, dir.path, prefix, reason = "manual") {
    dir.path <- path.expand(trimws(as.character(dir.path)))
    if (!nzchar(dir.path)) dir.path <- tempdir()
    dir.create(dir.path, recursive = TRUE, showWarnings = FALSE)
    idx <- sort(unique(as.integer(idx)))
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= n]
    pref <- sanitize.prefix(prefix)
    reason <- sanitize.prefix(reason)
    tag <- stamp.now()
    stem <- paste(pref, reason, tag, sep = "_")
    csv.file <- file.path(dir.path, paste0(stem, ".csv"))
    rds.file <- file.path(dir.path, paste0(stem, ".rds"))

    out.tbl <- data.frame(
      index = idx,
      rowname = if (!is.null(rownames(X))) as.character(rownames(X)[idx]) else NA_character_,
      stringsAsFactors = FALSE
    )
    out.tbl[[colnames(X)[1]]] <- X[idx, 1]
    out.tbl[[colnames(X)[2]]] <- X[idx, 2]
    out.tbl[[colnames(X)[3]]] <- X[idx, 3]
    write.csv(out.tbl, csv.file, row.names = FALSE)
    saveRDS(
      list(
        idx = idx,
        coords = X[idx, , drop = FALSE],
        plot.type = plot.type,
        generated.at = as.character(Sys.time())
      ),
      file = rds.file
    )
    list(
      csv = normalizePath(csv.file, mustWork = FALSE),
      rds = normalizePath(rds.file, mustWork = FALSE),
      dir = normalizePath(dir.path, mustWork = FALSE),
      n.selected = length(idx)
    )
  }

  hover.text <- function(idx) {
    rn <- if (!is.null(rownames(X))) rownames(X)[idx] else as.character(idx)
    paste0(
      "index=", idx,
      "<br>row=", rn,
      "<br>", colnames(X)[1], "=", signif(X[idx, 1], 5),
      "<br>", colnames(X)[2], "=", signif(X[idx, 2], 5),
      "<br>", colnames(X)[3], "=", signif(X[idx, 3], 5)
    )
  }

  build.plot <- function(selected.idx = integer(0)) {
    p <- plotly::plot_ly(source = source.id)

    if (length(idx.gray) > 0L) {
      p <- plotly::add_markers(
        p,
        x = X[idx.gray, 1],
        y = X[idx.gray, 2],
        z = X[idx.gray, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = point.size, color = "gray82", opacity = 0.25),
        customdata = idx.gray,
        text = hover.text(idx.gray),
        hoverinfo = "text",
        name = "gray",
        showlegend = FALSE
      )
    }

    if (length(idx.main) > 0L) {
      p <- plotly::add_markers(
        p,
        x = X[idx.main, 1],
        y = X[idx.main, 2],
        z = X[idx.main, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = point.size, color = colors[idx.main], opacity = 0.9),
        customdata = idx.main,
        text = hover.text(idx.main),
        hoverinfo = "text",
        name = "points",
        showlegend = FALSE
      )
    }

    if (length(selected.idx) > 0L) {
      selected.idx <- sort(unique(as.integer(selected.idx)))
      p <- plotly::add_markers(
        p,
        x = X[selected.idx, 1],
        y = X[selected.idx, 2],
        z = X[selected.idx, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = point.size + 3, color = "red", opacity = 1),
        customdata = selected.idx,
        text = hover.text(selected.idx),
        hoverinfo = "text",
        name = "selected",
        showlegend = FALSE
      )
    }

    if (length(end.labels) > 0L) {
      e.idx <- as.integer(names(end.labels))
      p <- plotly::add_trace(
        p,
        x = X[e.idx, 1],
        y = X[e.idx, 2],
        z = X[e.idx, 3],
        type = "scatter3d",
        mode = "text",
        text = unname(end.labels),
        textposition = "top center",
        hoverinfo = "skip",
        showlegend = FALSE,
        inherit = FALSE
      )
    }

    p <- plotly::layout(
      p,
      scene = list(
        xaxis = list(title = colnames(X)[1]),
        yaxis = list(title = colnames(X)[2]),
        zaxis = list(title = colnames(X)[3]),
        dragmode = "orbit"
      ),
      margin = list(l = 0, r = 0, b = 0, t = 0)
    )
    plotly::event_register(p, "plotly_click")
  }

  app <- shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::tags$h3(title),
      shiny::tags$div(
        style = "margin-bottom:8px;color:#444;",
        "Click points to toggle selection. Drag to rotate. Scroll to zoom."
      ),
      plotly::plotlyOutput("plot3d", height = "800px"),
      shiny::fluidRow(
        shiny::column(
          width = 4,
          shiny::actionButton("undo", "Undo"),
          shiny::tags$span(" "),
          shiny::actionButton("clear", "Clear"),
          shiny::tags$span(" "),
          shiny::actionButton("save_now", "Save Now"),
          shiny::tags$span(" "),
          shiny::actionButton("done", "Done (Save + Close)"),
          shiny::tags$span(" "),
          shiny::actionButton("cancel", "Cancel")
        ),
        shiny::column(
          width = 8,
          shiny::verbatimTextOutput("status")
        )
      ),
      if (isTRUE(allow.edit.save.dir)) {
        shiny::tagList(
          shiny::fluidRow(
            shiny::column(
              width = 12,
              shiny::textInput(
                "save_dir",
                "Save directory",
                value = as.character(save.dir),
                width = "min(96vw, 1800px)"
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(
              width = 6,
              shiny::textInput("save_prefix", "File prefix", value = as.character(save.prefix))
            )
          )
        )
      } else {
        shiny::fluidRow(
          shiny::column(
            width = 12,
            shiny::tags$small(
              style = "color:#666;",
              paste0("Save directory: ", path.expand(save.dir), " | prefix: ", save.prefix)
            )
          )
        )
      },
      shiny::tableOutput("preview")
    ),
    server = function(input, output, session) {
      selected.idx <- shiny::reactiveVal(integer(0))
      history <- shiny::reactiveVal(list())
      last.saved <- shiny::reactiveVal(NULL)

      get.save.dir <- function() {
        if (isTRUE(allow.edit.save.dir)) {
          v <- input$save_dir
          if (is.null(v) || !nzchar(v)) v <- as.character(save.dir)
          return(v)
        }
        as.character(save.dir)
      }
      get.save.prefix <- function() {
        if (isTRUE(allow.edit.save.dir)) {
          v <- input$save_prefix
          if (is.null(v) || !nzchar(v)) v <- as.character(save.prefix)
          return(v)
        }
        as.character(save.prefix)
      }

      shiny::observeEvent(plotly::event_data("plotly_click", source = source.id), {
        ed <- plotly::event_data("plotly_click", source = source.id)
        if (is.null(ed) || NROW(ed) < 1L || is.null(ed$customdata)) return()
        idx <- suppressWarnings(as.integer(ed$customdata[[1L]]))
        if (!is.finite(idx) || idx < 1L || idx > n) return()

        cur <- selected.idx()
        history(c(history(), list(cur)))
        if (idx %in% cur) {
          cur <- setdiff(cur, idx)
        } else {
          cur <- sort(unique(c(cur, idx)))
        }
        selected.idx(cur)
      }, ignoreNULL = TRUE)

      shiny::observeEvent(input$undo, {
        h <- history()
        if (length(h) < 1L) return()
        selected.idx(h[[length(h)]])
        history(h[seq_len(length(h) - 1L)])
      })

      shiny::observeEvent(input$clear, {
        history(c(history(), list(selected.idx())))
        selected.idx(integer(0))
      })

      shiny::observeEvent(input$save_now, {
        idx <- sort(unique(as.integer(selected.idx())))
        info <- try(write.selection.files(
          idx = idx,
          dir.path = get.save.dir(),
          prefix = get.save.prefix(),
          reason = "manual"
        ), silent = TRUE)
        if (inherits(info, "try-error")) {
          shiny::showNotification(
            paste0("Save failed: ", as.character(info)),
            type = "error",
            duration = 6
          )
          return()
        }
        last.saved(info)
        shiny::showNotification(
          paste0("Saved ", info$n.selected, " points to: ", info$csv),
          type = "message",
          duration = 4
        )
      })

      output$plot3d <- plotly::renderPlotly({
        build.plot(selected.idx())
      })

      output$status <- shiny::renderText({
        cur <- selected.idx()
        msg <- paste0("Selected: ", length(cur), " point(s)")
        if (!is.null(legend.info)) {
          msg <- paste0(msg, "\nColor bins/groups: ", length(legend.info))
        }
        ls <- last.saved()
        if (!is.null(ls)) {
          msg <- paste0(
            msg,
            "\nLast saved CSV: ", ls$csv,
            "\nLast saved RDS: ", ls$rds
          )
        }
        msg
      })

      output$preview <- shiny::renderTable({
        cur <- selected.idx()
        if (length(cur) < 1L) return(NULL)
        show.idx <- head(cur, 20L)
        data.frame(
          index = as.integer(show.idx),
          rowname = if (!is.null(rownames(X))) as.character(rownames(X)[show.idx]) else NA_character_,
          stringsAsFactors = FALSE
        )
      }, striped = TRUE, hover = TRUE)

      shiny::observeEvent(input$done, {
        idx <- sort(unique(as.integer(selected.idx())))
        saved.files <- last.saved()
        if (isTRUE(save.on.done)) {
          info <- try(write.selection.files(
            idx = idx,
            dir.path = get.save.dir(),
            prefix = get.save.prefix(),
            reason = "done"
          ), silent = TRUE)
          if (inherits(info, "try-error")) {
            shiny::showNotification(
              paste0("Done save failed: ", as.character(info), " (fix path or click Cancel)"),
              type = "error",
              duration = 8
            )
            return()
          }
          saved.files <- info
          last.saved(info)
        }
        shiny::stopApp(list(
          idx = idx,
          coords = X[idx, , drop = FALSE],
          n.selected = length(idx),
          plot.type = plot.type,
          cancelled = FALSE,
          saved.files = saved.files
        ))
      })

      shiny::observeEvent(input$cancel, {
        shiny::stopApp(empty.result)
      })
    }
  )

  result <- shiny::runApp(app, launch.browser = isTRUE(launch.browser))
  if (is.null(result)) result <- empty.result
  if (is.null(result$idx)) result$idx <- integer(0)
  result$idx <- sort(unique(as.integer(result$idx)))
  if (is.null(result$coords)) result$coords <- X[result$idx, , drop = FALSE]
  if (is.null(result$n.selected)) result$n.selected <- as.integer(length(result$idx))
  if (is.null(result$plot.type)) result$plot.type <- plot.type
  if (is.null(result$cancelled)) result$cancelled <- FALSE
  if (is.null(result$saved.files)) result$saved.files <- NULL

  invisible(result)
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

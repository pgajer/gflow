#' Plots 3D Points or Spheres Without Axes
#'
#' This function creates a 3D plot of points or spheres without axes labels
#' using the rgl package. It automatically switches between spheres and points
#' based on whether a radius is specified.
#'
#' @param X       matrix/data.frame with exactly 3 columns
#' @param radius  numeric or NULL; if non-NULL, draw spheres with this radius
#' @param col     color for points/spheres
#' @param size    point size when radius is NULL
#' @param axes,xlab,ylab,zlab standard rgl axis/label args
#' @param ...     passed to rgl::plot3d()
#'
#' @return Invisibly returns the rgl object ids from plot3d()
#' @export
plot3D.plain <- function(X,
                         radius = NULL,
                         col = "gray",
                         size = 3,
                         axes = FALSE, xlab = "", ylab = "", zlab = "",
                         ...) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl' for 3D visualization. ",
         "Install with install.packages('rgl').", call. = FALSE)
  }

  ## Headless/CI-safe; harmless on desktops
  use_null <- (!interactive()) ||
      identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
      (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
  old_opt <- options(rgl.useNULL = use_null)
  on.exit(options(old_opt), add = TRUE)

  ## ---- open/clear device ----
  ## Check if an rgl device is already open
  if (rgl::cur3d() == 0) {
      ## No device open, create a new one
      if (use_null) {
          ## Null device for headless environments
          rgl::open3d()
          on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
      } else {
          ## Interactive: create a large square window
          ## Get screen dimensions
          screen_info <- try(rgl::par3d("windowRect"), silent = TRUE)

          ## Calculate square window size (use ~80% of screen height for safety)
          if (inherits(screen_info, "try-error")) {
              ## Fallback if we can't get screen info
              window_size <- 800
          } else {
              ## Estimate available screen space
              ## par3d("windowRect") returns current window, not screen size
              ## Use a reasonable maximum
              window_size <- min(1200, 800)  ## Conservative default
          }

          rgl::open3d(windowRect = c(50, 50, 50 + window_size, 50 + window_size))
      }
  } else {
      ## Device already open, use it (but don't close it on exit)
      rgl::set3d(rgl::cur3d())
  }

  rgl::clear3d()

  if (!is.matrix(X) && !is.data.frame(X)) stop("X must be a matrix or data frame")
  if (ncol(X) != 3) stop("X must have exactly 3 columns")

  X <- as.matrix(X)
  storage.mode(X) <- "double"

  ids <- if (!is.null(radius)) {
    rgl::plot3d(X,
                axes = axes, xlab = xlab, ylab = ylab, zlab = zlab,
                type = "s", radius = radius, col = col, ...)
  } else {
    rgl::plot3d(X,
                axes = axes, xlab = xlab, ylab = ylab, zlab = zlab,
                type = "p", size = size, col = col, ...)
  }

  invisible(ids)
}

#' Creates 3D Plot with Continuous Color Coding
#'
#' Creates a 3D plot of a set of points defined by \code{X} color-coded by the values of \code{y}.
#'
#' @param X A matrix or data frame with 3 columns representing 3D coordinates.
#' @param y A numeric vector of values to be used for color coding.
#' @param subset A logical vector indicating which points to color-code by \code{y} values
#'   (default is NULL, which colors all points). When provided, points where \code{subset = TRUE}
#'   are drawn as colored spheres based on \code{y}, while points where \code{subset = FALSE}
#'   are drawn in \code{non.highlight.color} as background. Use this to emphasize a subset
#'   of points while showing the rest for context.
#' @param non.highlight.type Display style for non-highlighted points (where \code{subset = FALSE}):
#'   either "sphere" or "point" (default "sphere").
#' @param non.highlight.color Color for non-highlighted points (where \code{subset = FALSE});
#'   default "gray".
#' @param point.size Point size when \code{non.highlight.type="point"} (default 3).
#' @param legend.title Legend title (default "").
#' @param legend.cex Legend text size (default 2).
#' @param legend.side Side of the plot for the legend (default 3).
#' @param legend.line Line position for the legend (default 0.5).
#' @param radius Sphere radius (default NULL, that is, auto-scale).
#' @param quantize.method Method for quantizing y values: "uniform" or "quantile" (default "uniform").
#' @param quantize.wins.p Winsorization parameter for "uniform" (default 0.02).
#' @param quantize.round Logical; round quantization endpoints (default FALSE).
#' @param quantize.dig.lab Digits for quantization labels (default 2).
#' @param start,end Hue start/end for palette (defaults 1/6 and 0).
#' @param n.levels Number of color levels (default 10).
#'
#' @return Invisibly returns a list with \item{y.cols}{}, \item{y.col.tbl}{}, \item{legend.labs}{}.
#'
#' @seealso \code{\link{plot3D.plain}}, \code{\link[rgl]{plot3d}}, \code{\link[rgl]{spheres3d}}
#' @export
plot3D.cont <- function(X,
                        y,
                        subset = NULL,
                        non.highlight.type = "sphere",
                        non.highlight.color = "gray",
                        point.size = 3,
                        legend.title = "",
                        legend.cex = 2,
                        legend.side = 3,
                        legend.line = 0.5,
                        radius = NULL,
                        quantize.method = "uniform",
                        quantize.wins.p = 0.02,
                        quantize.round = FALSE,
                        quantize.dig.lab = 2,
                        start = 1/6, end = 0, n.levels = 10) {

    ## rgl is optional; error here because this function's purpose is plotting
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("This function requires the optional package 'rgl' for 3D visualization. ",
             "Install with install.packages('rgl').", call. = FALSE)
    }

    ## Headless/CI-safe; harmless on desktops
    use_null <- (!interactive()) ||
        identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
        (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
    old_opt <- options(rgl.useNULL = use_null)
    on.exit(options(old_opt), add = TRUE)

    ## ---- open/clear device ----
    ## Check if an rgl device is already open
    if (rgl::cur3d() == 0) {
        ## No device open, create a new one
        if (use_null) {
            ## Null device for headless environments
            rgl::open3d()
            on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
        } else {
            ## Interactive: create a large square window
            screen_info <- try(rgl::par3d("windowRect"), silent = TRUE)

            ## Calculate square window size (use ~80% of screen height for safety)
            if (inherits(screen_info, "try-error")) {
                ## Fallback if we can't get screen info
                window_size <- 800
            } else {
                ## Estimate available screen space
                ## par3d("windowRect") returns current window, not screen size
                ## Use a reasonable maximum
                window_size <- min(1200, 800)  ## Conservative default
            }

            rgl::open3d(windowRect = c(50, 50, 50 + window_size, 50 + window_size))
        }

    } else {
        ## Device already open, use it (but don't close it on exit)
        rgl::set3d(rgl::cur3d())
    }

    rgl::clear3d()

    ## if (rgl::cur3d() == 0) {
    ##     ## No device open, create a new one
    ##     if (use_null) {
    ##         ## Null device for headless environments
    ##         rgl::open3d()
    ##         on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
    ##     } else {
    ##         ## Interactive: create a large square window
    ##         ## Get screen dimensions
    ##         screen_info <- try(rgl::par3d("windowRect"), silent = TRUE)

    ##         ## Calculate square window size (use ~80% of screen height for safety)
    ##         if (inherits(screen_info, "try-error")) {
    ##             ## Fallback if we can't get screen info
    ##             window_size <- 800
    ##         } else {
    ##             ## Estimate available screen space
    ##             ## par3d("windowRect") returns current window, not screen size
    ##             ## Use a reasonable maximum
    ##             window_size <- min(1200, 800)  ## Conservative default
    ##         }

    ##         rgl::open3d(windowRect = c(50, 50, 50 + window_size, 50 + window_size))
    ##     }
    ## } else {
    ##     ## Device already open, use it (but don't close it on exit)
    ##     rgl::set3d(rgl::cur3d())
    ## }

    ## Validate inputs
    if (!is.matrix(X) && !is.data.frame(X)) stop("X must be a matrix or data frame")
    if (ncol(X) != 3) stop("X must have exactly 3 columns")
    if (!is.numeric(y)) stop("y must be a numeric vector")
    if (nrow(X) != length(y)) stop("Number of rows in X must match the length of y")
    if (!is.null(subset) && length(subset) != length(y)) {
        stop("Length of subset must match the length of y")
    }
    if (!non.highlight.type %in% c("sphere", "point")) {
        stop("non.highlight.type must be either 'sphere' or 'point'")
    }
    if (!is.numeric(point.size) || point.size <= 0) {
        stop("point.size must be a positive numeric value")
    }

    X <- as.matrix(X)
    storage.mode(X) <- "double"

    ## Quantize continuous y -> categories & colors
    q <- quantize.cont.var(
        y,
        method  = quantize.method,
        wins.p  = quantize.wins.p,
        round   = quantize.round,
        dig.lab = quantize.dig.lab,
        start = start, end = end, n.levels = n.levels
    )
    y.cat     <- q$x.cat
    y.col.tbl <- q$x.col.tbl
    y.cols    <- y.col.tbl[y.cat]

    ## Prepare subset; treat NAs as FALSE
    if (is.null(subset)) {
        subset <- rep(TRUE, length(y))
    } else {
        subset <- as.logical(subset)
        subset[is.na(subset)] <- FALSE
    }

    ## Choose a radius if none provided (scale-aware default)
    radius_local <- radius
    if (is.null(radius_local)) {
        rng <- apply(X, 2, function(v) diff(range(v, na.rm = TRUE)))
        ## fallback to a small positive number if data are near-constant
        radius_local <- max(1e-8, 0.01 * mean(rng))
    }

    ## Show which points have NA colors
    if (any(is.na(y.cols))) {
        cat("WARNING: Some points have NA colors!\n")
        cat("Indices with NA colors:", which(is.na(y.cols)), "\n")
    }

    ## Base layer: non-highlighted points
    if (any(!subset)) {
        if (identical(non.highlight.type, "sphere")) {
            plot3D.plain(
                X[!subset, , drop = FALSE],
                col      = non.highlight.color,
                radius   = radius_local,
                open_new = FALSE    # plot into the current device
            )
        } else {
            plot3D.plain(
                X[!subset, , drop = FALSE],
                col      = non.highlight.color,
                size     = point.size,
                open_new = FALSE
            )
        }
    } else {
        ## If no non-highlighted points, create an empty base so we can overlay
        plot3D.plain(
            X[FALSE, , drop = FALSE],
            col      = non.highlight.color,
            size     = point.size,
            open_new = FALSE
        )
    }

    ## Overlay: highlighted points (always spheres with data-driven colors)
    if (any(subset)) {
        rgl::spheres3d(X[subset, , drop = FALSE], col = y.cols[subset], radius = radius_local)
    }

    ## Legend: ensure all bins appear, including empty ones
    y.cat.freq <- table(y.cat)
    ## Fill missing bins with zero
    for (nm in names(y.col.tbl)) {
        if (!nm %in% names(y.cat.freq)) y.cat.freq[nm] <- 0L
    }
    ## Stable order by y.col.tbl names
    y.cat.freq <- y.cat.freq[names(y.col.tbl)]

    maxlen <- max(nchar(names(y.col.tbl)), 5) + 2
    legend.labs <- vapply(names(y.col.tbl), function(x) {
        sprintf(sprintf("%%-%ds%%5s", maxlen), x, sprintf("(%s)", y.cat.freq[[x]]))
    }, character(1))

    rgl::legend3d("topleft",
                  legend = legend.labs,
                  fill   = unname(y.col.tbl),
                  inset  = 0.05,
                  cex    = legend.cex,
                  title  = legend.title)

    invisible(list(
        y.cols     = y.cols,
        y.col.tbl  = y.col.tbl,
        legend.labs = legend.labs
    ))
}

#' Plot Output from lcor.1D()
#'
#' Creates a visualization of local correlation results from the lcor.1D() function
#'
#' @param r An output from lcor.1D() function.
#' @param with.CrI A logical parameter. If TRUE, credible intervals are plotted.
#' @param use.smoothed.lcors A logical parameter. If TRUE, smoothed local correlations will be used
#'   (if available) instead of raw local correlations.
#' @param title The title of the plot.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param lcor.col The color of local correlation line.
#' @param CrI.as.polygon A logical value parameter. If TRUE, credible interval region is drawn
#'   as a polygon. Otherwise, credible intervals are drawn as contour lines.
#' @param CrI.polygon.col The color of the credible interval polygon region.
#' @param CrI.line.col The color of the credible interval upper and lower limit lines.
#' @param add.box Set to TRUE if box() is to be called after plot(); useful when axes is FALSE.
#' @param add.x.axis Set to TRUE if axis(1) is to be called after plot().
#' @param add.y.axis Set to TRUE if axis(2) is to be called after plot().
#' @param ... Graphics parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating a plot.
#'
#' @details
#' This function visualizes the output from lcor.1D(), which typically contains local correlation
#' estimates along a one-dimensional variable. The function can display either raw or smoothed
#' local correlations, with optional credible intervals shown either as a shaded polygon or
#' as contour lines.
#'
#' @examples
#' \dontrun{
#' # Assuming 'result' is output from lcor.1D()
#' plotlcor.1D(result, with.CrI = TRUE, title = "Local Correlations")
#'
#' # Plot with credible intervals as lines instead of polygon
#' plotlcor.1D(result, with.CrI = TRUE, CrI.as.polygon = FALSE)
#' }
#'
#' @export
plotlcor.1D <- function(r, with.CrI = TRUE, use.smoothed.lcors = TRUE,
                       title = "", xlab = "", ylab = "Local Pearson correlation",
                       lcor.col = "blue", CrI.as.polygon = TRUE,
                       CrI.polygon.col = "gray95", CrI.line.col = "gray",
                       add.box = FALSE, add.x.axis = FALSE, add.y.axis = FALSE, ...) {

    if (use.smoothed.lcors && with.CrI) {
        plot(r$xg, r$smooth.lwcor.grid, type = "l", ylim = c(-1, 1), las = 1, col = lcor.col,
             ylab = ylab, xlab = xlab, main = title, ...)

        if (add.box) box()
        if (add.x.axis) axis(1)
        if (add.y.axis) axis(2, las = 1)

        if (CrI.as.polygon) {
            polygon(c(r$xg, rev(r$xg)),
                   c(r$smooth.lwcor.grid.CrI[1,], rev(r$smooth.lwcor.grid.CrI[2,])),
                   col = CrI.polygon.col, border = NA)
            lines(r$xg, r$smooth.lwcor.grid, col = lcor.col)
        } else {
            matlines(r$xg, t(r$smooth.lwcor.grid.CrI), col = CrI.line.col, lty = 1)
        }
    } else if (use.smoothed.lcors && !with.CrI) {
        plot(r$xg, r$smooth.lwcor.grid, type = "l", ylim = c(-1, 1), las = 1,
             col = lcor.col, ylab = ylab, xlab = xlab, main = title)
    } else if (with.CrI) {
        plot(r$xg, r$lwcor.grid, type = "l", ylim = c(-1, 1), las = 1,
             ylab = ylab, xlab = xlab, main = title)

        if (CrI.as.polygon) {
            polygon(c(r$xg, rev(r$xg)),
                   c(r$lwcor.grid.CrI[1,], rev(r$lwcor.grid.CrI[2,])),
                   col = CrI.polygon.col, border = NA)
            lines(r$xg, r$lwcor.grid, col = lcor.col)
        } else {
            matlines(r$xg, t(r$lwcor.grid.CrI), col = CrI.line.col, lty = 1)
        }
    } else {
        plot(r$xg, r$lwcor.grid, type = "l", ylim = c(-1, 1), las = 1,
             col = lcor.col, ylab = ylab, xlab = xlab, main = title)
    }

    abline(h = 0, col = "gray90")
    abline(h = -1, col = "gray")
    abline(h = 1, col = "gray")

    invisible(NULL)
}

#' Quantize Continuous Variable for Legend Display
#'
#' Quantizes a continuous variable and generates legend colors and labels for visualization
#'
#' @param y A numeric vector to be quantized.
#' @param quantize.method Method to use for quantization. Options are "uniform" (equal-width bins)
#'   or "quantile" (equal-frequency bins).
#' @param quantize.wins.p Winsorization parameter for the "uniform" method. Values beyond the
#'   p and 1-p quantiles are grouped together. Default is 0.02.
#' @param quantize.round Logical. If TRUE, rounds the endpoints of the range.
#' @param quantize.dig.lab Number of digits for labels in the cut() function.
#' @param start Start value for the color range in HSV color space (0-1).
#' @param end End value for the color range in HSV color space (0-1).
#' @param n.levels Number of discrete levels for quantization.
#'
#' @return A list containing:
#' \item{y.cat}{Factor of quantized values}
#' \item{y.col.tbl}{Named vector mapping categories to colors}
#' \item{y.cols}{Vector of colors for each input value}
#' \item{legend.labs}{Formatted labels for legend display}
#'
#' @details
#' This function is primarily used internally by plotting functions to convert continuous
#' variables into discrete categories with associated colors for visualization. The
#' winsorization parameter helps handle outliers in the uniform method.
#'
#' @examples
#' \dontrun{
#' y <- rnorm(100)
#' q <- quantize.for.legend(y, n.levels = 5)
#' plot(1:100, y, col = q$y.cols, pch = 19)
#' legend("topright", legend = q$legend.labs, fill = q$y.col.tbl)
#' }
#'
#' @export
quantize.for.legend <- function(y,
                               quantize.method = "uniform",
                               quantize.wins.p = 0.02,
                               quantize.round = FALSE,
                               quantize.dig.lab = 2,
                               start = 1/6, end = 0, n.levels = 10) {

    q <- quantize.cont.var(y,
                          method = quantize.method,
                          wins.p = quantize.wins.p,
                          round = quantize.round,
                          dig.lab = quantize.dig.lab,
                          start = start, end = end, n.levels = n.levels)

    y.cat <- q$x.cat
    y.col.tbl <- q$x.col.tbl
    y.cols <- y.col.tbl[y.cat]
    y.cat.freq <- table(y.cat)

    legend.labs <- sapply(names(y.col.tbl), FUN = function(x) {
        s <- sprintf("(%s)", y.cat.freq[x])
        sprintf("%-5s%5s", x, s)
    })

    list(y.cat = y.cat,
         y.col.tbl = y.col.tbl,
         y.cols = y.cols,
         legend.labs = legend.labs)
}

#' Create 2D Plot with Continuous Color Coding
#'
#' Creates a 2D plot of points color coded by the values of a continuous variable
#'
#' @param X A matrix or data.frame with at least 2 columns. Only the first 2 columns are used.
#' @param y A numeric vector of values for color coding.
#' @param legend.title Title for the legend.
#' @param legend.cex Character expansion factor for legend text.
#' @param legend.loc Location of the legend. If NULL, defaults to "topleft".
#' @param legend.line Line parameter for the legend.
#' @param legend.bty Box type for the legend.
#' @param use.layout Logical. If TRUE, uses layout() to position the legend separately.
#' @param layout.widths Numeric vector of length 2 specifying relative widths when use.layout is TRUE.
#' @param quantize.method Method for quantizing y values: "uniform" or "quantile".
#' @param quantize.wins.p Winsorization parameter for the "uniform" method.
#' @param quantize.round Logical. Whether to round the quantization endpoints.
#' @param quantize.dig.lab Number of digits for quantization labels.
#' @param start Start value for color range.
#' @param end End value for color range.
#' @param n.levels Number of color levels.
#' @param axes Logical. Whether to draw axes.
#' @param pch Plotting character.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns a list containing:
#' \item{y.cols}{Vector of colors assigned to each data point}
#' \item{y.col.tbl}{Color table mapping categories to colors}
#' \item{legend.labs}{Labels for the legend}
#' \item{legend.out}{Output from legend() function}
#'
#' @details
#' This function creates a 2D scatter plot where points are colored according to a continuous
#' variable. The continuous variable is discretized into categories for visualization.
#' The legend can be positioned using standard graphics layout or a custom layout approach.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(200), ncol = 2)
#' y <- runif(100)
#'
#' # Basic plot
#' plot2D.cont(X, y)
#'
#' # Plot with custom legend position and no axes
#' plot2D.cont(X, y, legend.title = "Values", axes = FALSE)
#' }
#'
#' @importFrom graphics plot.new grconvertX plot legend par layout
#' @export
plot2D.cont <- function(X,
                        y,
                        legend.title = "",
                        legend.cex = 0.7,
                        legend.loc = NULL,
                        legend.line = 0.5,
                        legend.bty = 'o',
                        use.layout = TRUE,
                        layout.widths = c(7, 1.5),
                        quantize.method = "uniform",
                        quantize.wins.p = 0.02,
                        quantize.round = FALSE,
                        quantize.dig.lab = 2,
                        start = 1/6, end = 0, n.levels = 10,
                        axes = FALSE,
                        pch = 20,
                        ...) {

    # Input validation
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame")
    }
    if (ncol(X) < 2) {
        stop("X must have at least 2 columns")
    }
    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }
    if (nrow(X) != length(y)) {
        stop("Number of rows in X must match the length of y")
    }

    q <- quantize.cont.var(y,
                          method = quantize.method,
                          wins.p = quantize.wins.p,
                          round = quantize.round,
                          dig.lab = quantize.dig.lab,
                          start = start, end = end, n.levels = n.levels)

    y.cat <- q$x.cat
    y.col.tbl <- q$x.col.tbl
    y.cols <- y.col.tbl[y.cat]
    y.cat.freq <- table(y.cat)

    legend.labs <- sapply(names(y.col.tbl), FUN = function(x) {
        s <- sprintf("(%s)", y.cat.freq[x])
        sprintf("%-5s%5s", x, s)
    })

    if (is.null(legend.loc)) {
        legend.loc <- "topleft"
    }

    opar <- par(no.readonly = TRUE)
    on.exit(par(opar), add = TRUE, after = FALSE)

    if (use.layout) {
        par(mar = c(0, 0, 0, 0))
        layout(matrix(c(1, 2), ncol = 2), widths = layout.widths)
        plot(X[, 1], X[, 2], col = y.cols, xlab = "", ylab = "", axes = axes, pch = pch, ...)

        # Setup for no margins on the legend
        par(mar = c(0, 0.3, 0, 0))
        plot.new()
        l <- legend("topleft", xpd = NA, legend = legend.labs, fill = y.col.tbl,
                   inset = 0.01, title = legend.title, cex = legend.cex)
    } else {
        plot.new()
        l <- legend(legend.loc, legend = legend.labs, fill = y.col.tbl,
                   title = legend.title, cex = legend.cex, bty = 'n', plot = FALSE)

        # Calculate right margin width in ndc
        w <- grconvertX(l$rect$w, to = 'ndc') - grconvertX(0, to = 'ndc')
        par(omd = c(0, 1 - w, 0, 1))
        plot(X[, 1], X[, 2], col = y.cols, xlab = "", ylab = "", axes = axes, pch = pch, ...)
        l <- legend(par('usr')[2], par('usr')[4], bty = legend.bty, xpd = NA,
                   legend = legend.labs, fill = y.col.tbl, inset = 0.05,
                   title = legend.title, cex = legend.cex)
    }

    invisible(list(y.cols = y.cols,
                   y.col.tbl = y.col.tbl,
                   legend.labs = legend.labs,
                   legend.out = l))
}

#' Display Cluster Labels in 3D Space (Headless/CRAN-safe)
#'
#' Shows cluster labels at the median center of each cluster in 3D space.
#' Computes centers regardless of plotting availability; plotting is optional.
#'
#' @param X A matrix or data.frame with 3 columns representing 3D coordinates.
#' @param cltr A vector of cluster assignments for the rows of \code{X}.
#'             \code{NA} entries are ignored.
#' @param cex Character expansion factor for cluster labels.
#' @param adj Adjustment values for label positioning (passed to \code{rgl::text3d}).
#' @param show.plot Logical; if \code{TRUE}, add labels to the current rgl scene
#'   (off-screen in headless CI). If rgl is not installed, plotting is skipped
#'   with a warning and centers are still returned.
#' @param open_new Logical; if \code{TRUE}, open a new rgl device for the labels
#'   and close it on exit. Defaults to \code{FALSE} so callers can layer on an
#'   existing scene.
#'
#' @return Invisibly returns a numeric matrix of size \eqn{K \times 3} containing
#'   per-cluster centers (component-wise medians when cluster size \eqn{\ge 2},
#'   or the lone point when size \eqn{=1}). Row names are the cluster labels.
#'
#' @details
#' This function is compute-first and plotting-optional. It is safe on CRAN/rhub
#' and in headless environments. When \code{show.plot = TRUE}, it renders
#' off-screen by setting \code{options(rgl.useNULL = TRUE)} locally. It never
#' requires \pkg{rgl} unless plotting is requested.
#'
#' @examples
#' \dontrun{
#' # Minimal example; runs even without rgl (no plotting).
#' set.seed(1)
#' X <- matrix(rnorm(300), ncol = 3)
#' cl <- sample(1:5, nrow(X), replace = TRUE)
#' centers <- show.cltrs(X, cl, show.plot = FALSE)
#'
#' # If rgl is available, render off-screen safely:
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
#'   # Suppose you've already drawn points via your helper:
#'   # plot3D.plain(X)
#'   show.cltrs(X, cl, cex = 1.1, adj = c(0.5, 1), show.plot = TRUE)
#' }
#' }
#'
#' @export
show.cltrs <- function(X, cltr, cex = 1, adj = c(0.5, 1),
                       show.plot = TRUE, open_new = FALSE) {
    ## --- Input checks & coercions (CRAN-safe) ---
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame", call. = FALSE)
    }
    X <- as.matrix(X)
    if (ncol(X) != 3L) {
        stop("X must have exactly 3 columns", call. = FALSE)
    }
    storage.mode(X) <- "double"

    if (length(cltr) != nrow(X)) {
        stop("Length of 'cltr' must match the number of rows in X", call. = FALSE)
    }

    ## Drop rows with NA cluster labels
    keep <- !is.na(cltr)
    if (!all(keep)) {
        X <- X[keep, , drop = FALSE]
        cltr <- cltr[keep]
    }
    if (nrow(X) == 0L) {
        warning("No non-NA cluster assignments; nothing to do.")
        return(invisible(matrix(numeric(0), nrow = 0L, ncol = 3,
                                dimnames = list(NULL, c("x", "y", "z")))))
    }

    ## Use character labels for rownames; preserve your freq-based ordering
    cltr_chr <- as.character(cltr)
    cltr_labs <- unique(names(sort(table(cltr_chr))))
    nCl <- length(cltr_labs)

    cltr_centers <- matrix(NA_real_, nrow = nCl, ncol = 3L,
                           dimnames = list(cltr_labs, c("x", "y", "z")))

    ## --- Compute per-cluster centers (medians for size >= 2, point itself for size == 1) ---
    for (j in seq_len(nCl)) {
        idx <- cltr_chr == cltr_labs[j]
        nj <- sum(idx)
        if (nj == 0L) next
        if (nj == 1L) {
            cltr_centers[j, ] <- X[which(idx)[1L], ]
        } else {
            ## median is robust; NA-safe because X has no NAs unless user provided them
            cltr_centers[j, ] <- apply(X[idx, , drop = FALSE], 2L, median, na.rm = TRUE)
        }
    }

    ## --- Optional plotting: headless-safe & gated on rgl availability ---
    if (isTRUE(show.plot)) {
        if (!requireNamespace("rgl", quietly = TRUE)) {
            warning("Package 'rgl' is not installed; skipping 3D visualization.")
        } else {

            use_null <- (!interactive()) ||
                identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
                (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
            old_opt <- options(rgl.useNULL = use_null)
            on.exit(options(old_opt), add = TRUE)

            if (isTRUE(open_new)) {
                rgl::open3d()
                if (use_null) {
                    ## Only close the device automatically if using null device
                    on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
                }
                rgl::clear3d()
            }
                                        # Add labels to the current scene (or the fresh one if open_new = TRUE)
            rgl::text3d(x = cltr_centers[, 1],
                        y = cltr_centers[, 2],
                        z = cltr_centers[, 3],
                        texts = rownames(cltr_centers),
                        font = 2L, cex = cex, adj = adj)
        }
    }

    invisible(cltr_centers)
}

#' Create 3D Plot with Cluster Visualization
#'
#' Creates a 3D plot of points with cluster assignments shown by different colors
#'
#' @param X A 3D matrix or data.frame with exactly 3 columns.
#' @param cltr A vector of cluster IDs (numeric or character) of length nrow(X).
#' @param cltr.col.tbl A named vector mapping cluster IDs to colors. If NULL, colors are assigned automatically.
#' @param ref.cltr A reference cluster ID (usually '0') that can be colored differently.
#' @param ref.cltr.color The color for the reference cluster.
#' @param show.ref.cltr Logical. Whether to display the reference cluster.
#' @param show.cltr.labels Logical. Whether to show cluster labels in the plot.
#' @param add Logical. Whether to add to an existing plot.
#' @param title Title of the plot.
#' @param cex.labs Size scaling parameter for cluster labels.
#' @param pal.type Palette type: "numeric", "brewer", or "mclust".
#' @param brewer.pal.n Number of colors in the RColorBrewer palette.
#' @param brewer.pal Name of the RColorBrewer palette.
#' @param show.legend Logical. Whether to show the legend.
#' @param sort.legend.labs.by.freq Logical. Sort legend labels by frequency.
#' @param sort.legend.labs.by.name Logical. Sort legend labels alphabetically.
#' @param filter.out.freq.0.cltrs Logical. Remove clusters with zero frequency from legend.
#' @param legend.title Title for the legend.
#' @param radius Numeric. Size of spheres at data points.
#' @param axes Logical. Whether to show axes.
#' @param xlab,ylab,zlab Axis labels.
#' @param ... Additional arguments passed to \code{\link[rgl]{plot3d}}.
#'
#' @return Invisibly returns a list containing:
#' \item{ids}{RGL object IDs}
#' \item{cltr.col.tbl}{Color table used for clusters}
#' \item{cltr.labs}{Cluster labels}
#' \item{legend.cltr.labs}{Legend labels}
#' \item{cltr.centers}{Matrix of cluster centers (if show.cltr.labels = TRUE)}
#'
#' @details
#' This function provides comprehensive cluster visualization in 3D space with automatic
#' color assignment, legend generation, and various customization options. It supports
#' highlighting specific clusters and different color palette options.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
#'   set.seed(123)
#'   X <- matrix(rnorm(300), ncol = 3)
#'   cltr <- sample(c("A", "B", "C"), 100, replace = TRUE)
#'   plot3D.cltrs(X, cltr, show.cltr.labels = TRUE)
#' }
#' }
#' @importFrom grDevices hcl.colors rainbow
#' @export
#' @export
plot3D.cltrs <- function(X,
                         cltr = NULL,
                         cltr.col.tbl = NULL,
                         ref.cltr = NULL,
                         ref.cltr.color = 'gray',
                         show.ref.cltr = TRUE,
                         show.cltr.labels = FALSE,
                         add = FALSE,
                         title = "",
                         cex.labs = 2,
                         pal.type = "numeric",
                         brewer.pal.n = 3,
                         brewer.pal = "Set1",
                         show.legend = TRUE,
                         sort.legend.labs.by.freq = FALSE,
                         sort.legend.labs.by.name = FALSE,
                         filter.out.freq.0.cltrs = TRUE,
                         legend.title = NULL,
                         radius = NA,
                         axes = FALSE,
                         xlab = "",
                         ylab = "",
                         zlab = "",
                         ...) {

    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("This function requires the optional package 'rgl' for 3D visualization. ",
             "Install with install.packages('rgl').", call. = FALSE)
    }

    ## Headless/CI-safe; harmless on desktops
    use_null <- (!interactive()) ||
        identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
        (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
    old_opt <- options(rgl.useNULL = use_null)
    on.exit(options(old_opt), add = TRUE)

    ## ---- open/clear device ----
    ## Check if an rgl device is already open
    if (rgl::cur3d() == 0) {
        ## No device open, create a new one
        if (use_null) {
            ## Null device for headless environments
            rgl::open3d()
            on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
        } else {
            ## Interactive: create a large square window
            ## Get screen dimensions
            screen_info <- try(rgl::par3d("windowRect"), silent = TRUE)

            ## Calculate square window size (use ~80% of screen height for safety)
            if (inherits(screen_info, "try-error")) {
                ## Fallback if we can't get screen info
                window_size <- 800
            } else {
                ## Estimate available screen space
                ## par3d("windowRect") returns current window, not screen size
                ## Use a reasonable maximum
                window_size <- min(1200, 800)  ## Conservative default
            }

            rgl::open3d(windowRect = c(50, 50, 50 + window_size, 50 + window_size))
        }
    } else {
        ## Device already open, use it (but don't close it on exit)
        rgl::set3d(rgl::cur3d())
    }

    rgl::clear3d()

    if (!is.matrix(X) && !is.data.frame(X)) stop("X must be a matrix or data frame")
    if (ncol(X) != 3) stop("X must have exactly 3 columns")
    if (!is.null(cltr) && length(cltr) != nrow(X)) stop("Length of cltr must match nrow(X)")

    legend.cltr.labs <- NULL
    ids <- NULL
    cltr.centers <- NULL

    x <- X[, 1]; y <- X[, 2]; z <- X[, 3]

    if (!is.null(cltr)) {
        cltr <- as.character(cltr)
        cltr.labs <- names(sort(table(cltr)))
        nCl <- length(cltr.labs)

        if (is.null(cltr.col.tbl)) {
            if (nCl <= 8) pal.type <- "numeric"
            if (nCl > 19) {
                brewer.pal <- "Spectral"; brewer.pal.n <- 11
                cltr.labs <- sample(cltr.labs)
            }

            if (pal.type == "brewer") {
                if (requireNamespace("RColorBrewer", quietly = TRUE)) {
                    col.pal <- grDevices::colorRampPalette(
                                              rev(RColorBrewer::brewer.pal(brewer.pal.n, brewer.pal))
                                          )
                    cltr.col.tbl <- col.pal(nCl)
                } else {
                    warning("RColorBrewer not installed; using hcl.colors('Spectral') fallback.")
                    cltr.col.tbl <- grDevices::hcl.colors(nCl, palette = "Spectral")
                }
            } else if (pal.type == "numeric") {
                cltr.col.tbl <- if (nCl < 8) 2:(nCl + 1) else seq_len(nCl)
            } else if (pal.type == "mclust") {
                if (!requireNamespace("mclust", quietly = TRUE)) {
                    warning("mclust not installed; using hcl.colors('Spectral') fallback.")
                    cltr.col.tbl <- grDevices::hcl.colors(nCl, palette = "Spectral")
                } else {
                    cltr.col.tbl <- mclust::mclust.options("classPlotColors")[seq_len(nCl)]
                }
            } else {
                stop(sprintf("Unrecognized pal.type='%s'", pal.type))
            }
            names(cltr.col.tbl) <- cltr.labs
        }

        if (!is.null(ref.cltr)) {
            ref.cltr <- as.character(ref.cltr)
            if (!(ref.cltr %in% cltr.labs)) warning(ref.cltr, " is not in cltr.labs")
            cltr.col.tbl[ref.cltr] <- ref.cltr.color
        }

        cltr.cols <- cltr.col.tbl[cltr]
        if (nCl == 1 && !add) {
            ids <- rgl::plot3d(x, y, z, axes = axes, xlab = xlab, ylab = ylab, zlab = zlab, main = title, ...)
        } else {
            if (!is.na(radius)) {
                if (!add) {
                    if (is.null(ref.cltr)) {
                        ids <- rgl::plot3d(x, y, z, col = cltr.cols, axes = axes, xlab = xlab, ylab = ylab, zlab = zlab,
                                           main = title, type = "s", radius = radius, ...)
                    } else {
                        idx <- cltr != ref.cltr
                        ids <- rgl::plot3d(x[idx], y[idx], z[idx], col = cltr.cols[idx], axes = axes,
                                           xlab = xlab, ylab = ylab, zlab = zlab, main = title, type = "s", radius = radius, ...)
                        if (show.ref.cltr) rgl::spheres3d(x[!idx], y[!idx], z[!idx], col = ref.cltr.color, radius = radius)
                    }
                }
            } else if (!add) {
                ids <- rgl::plot3d(x, y, z, col = cltr.cols, axes = axes, xlab = xlab, ylab = ylab, zlab = zlab, main = title, ...)
            }

            if (show.cltr.labels) {
                cltr.centers <- matrix(nrow = length(cltr.labs), ncol = 3,
                                       dimnames = list(cltr.labs, c("x","y","z")))
                for (j in seq_along(cltr.labs)) {
                    idx <- cltr == cltr.labs[j]
                    idx[is.na(idx)] <- FALSE
                    s <- sum(idx)
                    if (s > 1) cltr.centers[j, ] <- apply(X[idx, , drop = FALSE], 2, median)
                    else if (s == 1) cltr.centers[j, ] <- as.numeric(X[idx, ])
                }
                rgl::text3d(x = cltr.centers[,1], y = cltr.centers[,2], z = cltr.centers[,3],
                            texts = rownames(cltr.centers), font = 2, cex = cex.labs, adj = c(0.5, 1))
            }

            cltr.freq <- table(cltr)
            if (!is.null(ref.cltr)) cltr.freq <- cltr.freq[setdiff(names(cltr.freq), ref.cltr)]

            if (isTRUE(show.legend)) {
                for (xnm in names(cltr.col.tbl)) if (!(xnm %in% names(cltr.freq))) cltr.freq[xnm] <- 0L
                if (isTRUE(sort.legend.labs.by.freq)) {
                    o <- order(cltr.freq[names(cltr.col.tbl)], decreasing = TRUE); cltr.col.tbl <- cltr.col.tbl[o]
                } else if (isTRUE(sort.legend.labs.by.name)) {
                    o <- order(names(cltr.col.tbl)); cltr.col.tbl <- cltr.col.tbl[o]
                }
                if (isTRUE(filter.out.freq.0.cltrs)) {
                    cltr.freq <- cltr.freq[names(cltr.col.tbl)]
                    keep <- cltr.freq != 0L
                    cltr.col.tbl <- cltr.col.tbl[keep]; cltr.freq <- cltr.freq[keep]
                }
                maxlen <- max(nchar(names(cltr.col.tbl))) + 2
                legend.cltr.labs <- vapply(names(cltr.col.tbl), function(xnm) {
                    sprintf(sprintf("%%-%ds%%5s", maxlen), xnm, sprintf("(%s)", cltr.freq[xnm]))
                }, character(1))
                rgl::legend3d("topleft", legend = legend.cltr.labs, fill = cltr.col.tbl, inset = 0.05,
                              cex = 1.5, title = legend.title)
            } else {
                legend.cltr.labs <- NULL
                try(rgl::legend3d(), silent = TRUE)
            }
        }
    } else {
        ids <- rgl::plot3d(X, axes = axes, xlab = xlab, ylab = ylab, zlab = zlab, main = title, ...)
        cltr.col.tbl <- NULL; cltr.labs <- NULL
        try(rgl::legend3d(), silent = TRUE)
    }

    invisible(list(ids = ids,
                   cltr.col.tbl = cltr.col.tbl,
                   cltr.labs = if (!is.null(cltr)) cltr.labs else NULL,
                   legend.cltr.labs = legend.cltr.labs,
                   cltr.centers = cltr.centers))
}

#' Highlight a Specific Cluster in 3D Space
#'
#' Shows a specified cluster with emphasis while displaying other clusters in gray
#'
#' @param cl Cluster ID to highlight.
#' @param cltr Vector of cluster IDs of length nrow(X).
#' @param X A matrix or data.frame with 3 columns representing 3D coordinates.
#' @param cl.radius Radius of spheres for the highlighted cluster.
#' @param show.ref.cltr Logical. Whether to show non-highlighted clusters.
#' @param show.labels Logical. Whether to show row names of highlighted points.
#' @param adj Adjustment parameter for text positioning.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns NULL. The function is called for its side effect.
#'
#' @details
#' This function creates a 3D plot where one specific cluster is highlighted in red
#' while all other clusters are shown in gray. Optionally, labels can be shown for
#' the highlighted cluster points.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(300), ncol = 3)
#' cltr <- sample(1:5, 100, replace = TRUE)
#'
#' # Highlight cluster 3
#' show.3d.cl(3, cltr, X)
#'
#' # Show with labels
#' rownames(X) <- paste0("Point", 1:100)
#' show.3d.cl(3, cltr, X, show.labels = TRUE)
#' }
#'
#' @export
show.3d.cl <- function(cl, cltr, X, cl.radius = 0.0001, show.ref.cltr = TRUE,
                       show.labels = FALSE, adj = c(1.3, 0), ...) {

    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame")
    }
    if (ncol(X) != 3) {
        stop("X must have exactly 3 columns")
    }
    if (length(cltr) != nrow(X)) {
        stop("Length of cltr must match number of rows in X")
    }

    idx <- cltr == cl
    lab <- paste0("not ", cl)
    loc.cltr <- ifelse(idx, cl, lab)

    plot3D.cltrs(X, loc.cltr, ref.cltr = lab, show.ref.cltr = show.ref.cltr, ...)
    rgl::spheres3d(X[idx, , drop = FALSE], radius = cl.radius, col = 'red')

    if (show.labels && !is.null(rownames(X))) {
        rgl::text3d(X[idx, , drop = FALSE], texts = rownames(X)[idx], adj = adj)
    }

    invisible(NULL)
}

#' Plot a Specific Cluster with Custom Colors
#'
#' Creates a 3D plot highlighting one or more specific clusters
#'
#' @param cl Cluster ID(s) to highlight. Can be a single value or vector.
#' @param cltr Vector of cluster IDs.
#' @param X A matrix or data.frame with 3 columns representing 3D coordinates.
#' @param ... Additional arguments passed to \code{\link{plot3D.cltrs}}.
#'
#' @return Invisibly returns the output from plot3D.cltrs.
#'
#' @details
#' This function creates a 3D plot where specified clusters are highlighted with
#' distinct colors while other clusters are shown in gray shades. If multiple
#' clusters are specified, they are colored using the mclust color scheme.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(300), ncol = 3)
#' cltr <- sample(1:5, 100, replace = TRUE)
#'
#' # Highlight single cluster
#' plot3D.cl(2, cltr, X)
#'
#' # Highlight multiple clusters
#' plot3D.cl(c(2, 4), cltr, X)
#' }
#'
#' @export
plot3D.cl <- function(cl, cltr, X, ...) {

    cl <- as.character(cl)
    all <- as.character(sort(unique(cltr)))
    other <- setdiff(all, cl)

    cltr.col.tbl <- rep("gray", length(all))
    names(cltr.col.tbl) <- all

    if (length(cl) > 1) {
        if (!requireNamespace("mclust", quietly = TRUE)) {
            # Use default colors if mclust not available
            cltr.col.tbl[cl] <- 2:(length(cl) + 1)
        } else {
            cltr.col.tbl[cl] <- mclust::mclust.options("classPlotColors")[seq_along(cl)]
        }
    } else {
        cltr.col.tbl[cl] <- "red"
    }

    cltr.col.tbl[other] <- grDevices::gray.colors(length(other))

    plot3D.cltrs(X = X, cltr = cltr, cltr.col.tbl = cltr.col.tbl, ...)
}

#' Add 3D Line Segments for Binary Variable
#'
#' Adds vertical line segments to a 3D plot at positions where a binary variable equals 1
#'
#' @param X A matrix or data.frame with 3 columns representing 3D coordinates.
#' @param y A named binary (0/1) vector.
#' @param offset Numeric vector of length 3 specifying the offset for line segments.
#' @param with.labels Logical. Whether to show labels for y=1 positions.
#' @param lab.tbl Named vector mapping sample IDs to labels.
#' @param lab.adj Adjustment parameter for label positioning.
#' @param lab.cex Character expansion factor for labels.
#' @param C Scaling factor for label position relative to stick center.
#' @param ... Additional arguments passed to \code{\link[rgl]{segments3d}}.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function adds vertical line segments (sticks) to an existing 3D plot
#' at positions where the binary variable y equals 1. Optionally, labels can
#' be added above the sticks.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(300), ncol = 3)
#' y <- sample(0:1, 100, replace = TRUE)
#' names(y) <- rownames(X) <- paste0("Sample", 1:100)
#'
#' plot3D.plain(X)
#' bin.segments3d(X, y, offset = c(0, 0, 0.1))
#' }
#'
#' @export
bin.segments3d <- function(X, y, offset, with.labels = TRUE, lab.tbl = NULL,
                          lab.adj = c(0, 0), lab.cex = 1, C = 1, ...) {

    if (is.null(names(y))) {
        stop("y must be a named vector")
    }

    if (!all(rownames(X) == names(y))) {
        stop("rownames(X) must match names(y)")
    }

    if (!all(y %in% c(0, 1, NA))) {
        stop("y must be a binary (0/1) vector")
    }

    ids <- rownames(X)[y == 1]

    for (id in ids) {
        x <- X[id, ]
        M <- rbind(x - offset, x + offset)
        rgl::segments3d(M, ...)

        if (with.labels && !is.null(lab.tbl)) {
            rgl::text3d(x[1] + C * offset[1],
                       x[2] + C * offset[2],
                       x[3] + C * offset[3],
                       texts = lab.tbl[id],
                       adj = lab.adj,
                       cex = lab.cex)
        }
    }

    invisible(NULL)
}

#' Add 3D Line Segments Color-Coded by Continuous Variable
#'
#' Adds line segments to a 3D plot with colors determined by a continuous variable
#'
#' @param X A matrix or data.frame with 3 columns representing 3D coordinates.
#' @param y A named continuous numeric vector.
#' @param offset Numeric vector of length 3 specifying the offset for line segments.
#' @param ... Additional arguments passed to \code{\link[rgl]{segments3d}}.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function adds line segments to an existing 3D plot where each segment's
#' color is determined by the rank of the corresponding y value. Colors range
#' from red to violet following the rainbow spectrum.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(300), ncol = 3)
#' y <- runif(100)
#' names(y) <- rownames(X) <- paste0("Sample", 1:100)
#'
#' plot3D.plain(X)
#' cont.segments3d(X, y, offset = c(0, 0, 0.1))
#' }
#'
#' @importFrom grDevices rainbow
#' @export
cont.segments3d <- function(X, y, offset, ...) {

    if (!all(rownames(X) == names(y))) {
        stop("rownames(X) must match names(y)")
    }

    if (!is.numeric(y)) {
        stop("y must be numeric")
    }

    y.rk <- rank(y)
    col.tbl <- grDevices::rainbow(length(y), start = 1/6, end = 0)
    y.col <- col.tbl[y.rk]

    for (i in seq(nrow(X))) {
        x <- X[i, ]
        M <- rbind(x - offset, x + offset)
        rgl::segments3d(M, col = y.col[i], ...)
    }

    invisible(NULL)
}

#' Plot Median Absolute Error
#'
#' Creates a plot of Median Absolute Error (MAE) with error bars
#'
#' @param mae.mean Vector of MAE means.
#' @param mae.mad Vector of MAE median absolute deviations.
#' @param ylab Label for y-axis.
#' @param xlab Label for x-axis.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function creates a line plot with error bars showing the median absolute
#' error across different numbers of nearest neighbors. The error bars represent
#' +/-0.5 * MAD (median absolute deviation).
#'
#' @examples
#' \dontrun{
#' mae.mean <- c(NA, 0.5, 0.4, 0.35, 0.33, 0.32)
#' mae.mad <- c(NA, 0.1, 0.08, 0.07, 0.06, 0.06)
#' mae.plot(mae.mean, mae.mad)
#' }
#'
#' @seealso \code{\link{vert.error.bar}}
#' @export
mae.plot <- function(mae.mean, mae.mad,
                    ylab = "Median Absolute Error",
                    xlab = "Number of Nearest Neighbors") {

    if (length(mae.mean) != length(mae.mad)) {
        stop("mae.mean and mae.mad must have the same length")
    }

    r <- which(!is.na(mae.mean))
    if (length(r) == 0) {
        stop("mae.mean contains only NA values")
    }

    k.min <- min(r)
    k.max <- max(r)

    plot(k.min:k.max, mae.mean[k.min:k.max], type = 'b',
         xlab = xlab, ylab = ylab, las = 1,
         ylim = c(min(mae.mean - mae.mad/2, na.rm = TRUE),
                  max(mae.mean + mae.mad/2, na.rm = TRUE)))

    for (k in k.min:k.max) {
        if (!is.na(mae.mean[k]) && !is.na(mae.mad[k])) {
            vert.error.bar(k, mae.mean[k] - mae.mad[k]/2,
                          mae.mean[k] + mae.mad[k]/2, lwd = 1)
        }
    }

    invisible(NULL)
}

#' Map Sample IDs to 3D Space
#'
#' Highlights a subset of samples in 3D space with different colors and sizes
#'
#' @param S Vector of sample IDs to highlight.
#' @param X Matrix or data.frame with 3 columns representing 3D coordinates.
#' @param radius Radius of spheres for highlighted samples.
#' @param col Color for highlighted samples.
#' @param legend.title Title for the legend.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function creates a 3D plot where samples in set S are highlighted
#' with larger spheres in a specified color, while other samples are shown
#' in gray.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(300), ncol = 3)
#' rownames(X) <- paste0("Sample", 1:100)
#' S <- paste0("Sample", c(1, 5, 10, 15, 20))
#'
#' map.S.to.X(S, X, radius = 0.075, col = 'red')
#' }
#'
#' @export
map.S.to.X <- function(S, X, radius = 0.075, col = 'red', legend.title = NULL) {

    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame")
    }
    if (ncol(X) != 3) {
        stop("X must have exactly 3 columns")
    }

    cn <- intersect(S, rownames(X))
    if (length(cn) == 0) {
        warning("No samples from S found in rownames(X)")
    }

    ind <- numeric(nrow(X))
    names(ind) <- rownames(X)
    ind[cn] <- 1

    # Assuming bin.col.tbl is defined elsewhere or using default
    bin.col.tbl <- c("0" = "gray", "1" = col)

    if (is.null(legend.title)) {
        plot3D.cltrs(X, ind, cltr.col.tbl = bin.col.tbl)
    } else {
        plot3D.cltrs(X, ind, cltr.col.tbl = bin.col.tbl, legend.title = legend.title)
    }

    if (length(cn) > 0) {
        rgl::spheres3d(X[cn, , drop = FALSE], col = col, radius = radius)
    }

    invisible(NULL)
}

#' Add Minimal Spanning Tree to 3D Plot
#'
#' Adds edges of a minimal spanning tree to an existing 3D plot
#'
#' @param X Matrix with 3 columns representing 3D coordinates.
#' @param T.edges Matrix with 2 columns containing start and end indices of edges.
#' @param col Color of the edges.
#' @param lwd Line width of the edges.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(30), ncol = 3)
#' # Assuming T.edges is computed from a minimal spanning tree algorithm
#' T.edges <- matrix(c(1,2, 2,3, 3,4, 4,5), ncol = 2, byrow = TRUE)
#'
#' plot3D.plain(X)
#' plot3D.tree(X, T.edges)
#' }
#'
#' @export
plot3D.tree <- function(X, T.edges, col = "gray", lwd = 1) {

    if (!is.matrix(X) || ncol(X) != 3) {
        stop("X must be a matrix with 3 columns")
    }

    if (!is.matrix(T.edges) || ncol(T.edges) != 2) {
        stop("T.edges must be a matrix with 2 columns")
    }

    max_idx <- max(T.edges)
    if (max_idx > nrow(X)) {
        stop("Edge indices exceed number of rows in X")
    }

    for (i in seq(nrow(T.edges))) {
        s <- T.edges[i, 1]
        e <- T.edges[i, 2]
        rgl::segments3d(x = X[c(s, e), 1],
                       y = X[c(s, e), 2],
                       z = X[c(s, e), 3],
                       col = col, lwd = lwd)
    }

    invisible(NULL)
}

#' Add Minimal Spanning Tree to 2D Plot
#'
#' Adds edges of a minimal spanning tree to an existing 2D plot
#'
#' @param X Matrix with at least 2 columns. Only first 2 columns are used.
#' @param T.edges Matrix with 2 columns containing start and end indices of edges.
#' @param col Color of the edges.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(20), ncol = 2)
#' T.edges <- matrix(c(1,2, 2,3, 3,4, 4,5), ncol = 2, byrow = TRUE)
#'
#' plot(X, pch = 19)
#' plot2D.tree(X, T.edges)
#' }
#'
#' @export
plot2D.tree <- function(X, T.edges, col = "gray") {

    if (!is.matrix(X) || ncol(X) < 2) {
        stop("X must be a matrix with at least 2 columns")
    }

    if (!is.matrix(T.edges) || ncol(T.edges) != 2) {
        stop("T.edges must be a matrix with 2 columns")
    }

    max_idx <- max(T.edges)
    if (max_idx > nrow(X)) {
        stop("Edge indices exceed number of rows in X")
    }

    for (i in seq(nrow(T.edges))) {
        s <- T.edges[i, 1]
        e <- T.edges[i, 2]
        segments(X[s, 1], X[s, 2], X[e, 1], X[e, 2], col = col)
    }

    invisible(NULL)
}

#' Plot Path in 3D Graph
#'
#' Plots a sequence of edges forming a path in a graph with 3D vertex positions
#'
#' @param s Vector of vertex indices defining the path.
#' @param V Matrix of vertex positions with 3 columns.
#' @param edge.col Color of the path edges.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function draws line segments connecting consecutive vertices in the
#' specified path sequence.
#'
#' @examples
#' \dontrun{
#' V <- matrix(rnorm(15), ncol = 3)
#' s <- c(1, 3, 2, 5, 4)
#'
#' plot3D.plain(V)
#' plot3D.path(s, V, edge.col = "red")
#' }
#'
#' @export
plot3D.path <- function(s, V, edge.col = "gray") {

    if (!is.numeric(s) || length(s) < 2) {
        stop("s must be a numeric vector with at least 2 elements")
    }

    if (!is.matrix(V) || ncol(V) != 3) {
        stop("V must be a matrix with 3 columns")
    }

    if (max(s) > nrow(V)) {
        stop("Path indices exceed number of vertices")
    }

    for (i in 1:(length(s) - 1)) {
        M <- rbind(V[s[i], ], V[s[i + 1], ])
        rgl::segments3d(M, col = edge.col)
    }

    invisible(NULL)
}

#' Plot Geodesic Path in 3D
#'
#' Plots the shortest path between two vertices in a graph
#'
#' @param S.3d Matrix of 3D vertex positions.
#' @param i1 Index of the starting vertex.
#' @param i2 Index of the ending vertex.
#' @param G An igraph object representing the graph structure.
#' @param X Original state space (currently unused but kept for compatibility).
#' @param edge.col Color of the geodesic path.
#'
#' @return Invisibly returns the sequence of vertex indices in the shortest path.
#'
#' @details
#' This function finds the shortest path between two vertices using igraph's
#' shortest_paths function and visualizes it in 3D space.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' # Create a simple graph
#' G <- make_ring(10)
#' S.3d <- matrix(rnorm(30), ncol = 3)
#'
#' plot3D.plain(S.3d)
#' plot3D.geodesic(S.3d, 1, 5, G, S.3d)
#' }
#'
#' @importFrom igraph shortest_paths
#' @export
plot3D.geodesic <- function(S.3d, i1, i2, G, X, edge.col = "gray") {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required. Please install it.")
    }

    if (!is.matrix(S.3d) || ncol(S.3d) != 3) {
        stop("S.3d must be a matrix with 3 columns")
    }

    s <- igraph::shortest_paths(G, from = i1, to = i2)$vpath[[1]]

    if (length(s) == 0) {
        warning("No path found between vertices ", i1, " and ", i2)
        return(invisible(numeric(0)))
    }

    rgl::lines3d(S.3d[s, ], col = edge.col)

    invisible(as.numeric(s))
}

#' Animate Path Along Trajectory
#'
#' Function for animation of a smoothed trajectory
#'
#' @param i Index in the Epath matrix.
#' @param Epath Matrix of smoothed trajectory positions.
#' @param radius Radius of the animated sphere.
#' @param sphere.col Color of the animated sphere.
#' @param path.col Color of the path trail.
#' @param path.lwd Width of the path trail.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function is designed to be used in animation loops. It draws a sphere
#' at the current position and adds a line segment from the previous position.
#' Note: This function assumes 'sphere.id' exists in the calling environment.
#'
#' @examples
#' \dontrun{
#' # This function is typically used within an animation loop
#' Epath <- matrix(rnorm(30), ncol = 3)
#' sphere.id <- NULL
#' for (i in 2:nrow(Epath)) {
#'   path.fn(i, Epath)
#'   Sys.sleep(0.1)
#' }
#' }
#'
#' @export
path.fn <- function(i, Epath, radius = 0.0008, sphere.col = "red",
                    path.col = "blue", path.lwd = 5) {

    if (i < 2 || i > nrow(Epath)) {
        stop("i must be between 2 and nrow(Epath)")
    }

    rgl::par3d(skipRedraw = TRUE)
    on.exit(rgl::par3d(skipRedraw = FALSE))

    # Check if sphere.id exists in parent environment
    if (exists("sphere.id", envir = parent.frame()) &&
        !is.null(get("sphere.id", envir = parent.frame()))) {
        try(rgl::pop3d(id = get("sphere.id", envir = parent.frame())), silent = TRUE)
    }

    x1 <- Epath[i - 1, ]
    x2 <- Epath[i, ]
    M <- rbind(x1, x2)
    rgl::segments3d(M, col = path.col, lwd = path.lwd)
    sphere.id <- rgl::spheres3d(x1, col = sphere.col, radius = radius)

    # Store sphere.id in parent environment for next iteration
    assign("sphere.id", sphere.id, envir = parent.frame())

    invisible(NULL)
}

#' Add Gradient Arrows to 3D Plot
#'
#' Adds gradient vectors as 3D arrows to an existing plot
#'
#' @param ids Vector of row names to display gradients for.
#' @param S Matrix of 3D positions.
#' @param grad.ED Matrix of gradient vectors (same dimensions as S).
#' @param C Scaling constant for gradient visualization.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function adds arrows showing gradient directions at specified points.
#' The arrows are scaled by factor C for better visualization.
#'
#' @examples
#' \dontrun{
#' S <- matrix(rnorm(30), ncol = 3)
#' rownames(S) <- paste0("Point", 1:10)
#' grad.ED <- matrix(rnorm(30) * 0.1, ncol = 3)
#' rownames(grad.ED) <- rownames(S)
#'
#' plot3D.plain(S)
#' add.grad.ED.arrows(c("Point1", "Point5"), S, grad.ED, C = 2.5)
#' }
#'
#' @export
add.grad.ED.arrows <- function(ids, S, grad.ED, C = 2.5) {

    if (nrow(S) != nrow(grad.ED)) {
        stop("S and grad.ED must have the same number of rows")
    }

    if (ncol(S) != ncol(grad.ED)) {
        stop("S and grad.ED must have the same number of columns")
    }

    if (ncol(S) != 3) {
        stop("S must have exactly 3 columns")
    }

    comm.ids <- intersect(ids, rownames(grad.ED))
    if (length(comm.ids) < length(ids)) {
        warning("Some ids are not present in rownames(grad.ED)")
    }

    for (id in comm.ids) {
        if (sum(abs(grad.ED[id, ])) > 0) {
            v0 <- as.vector(S[id, ])
            v1 <- v0 + C * as.vector(grad.ED[id, ])
            rgl::arrow3d(v0, v1, add = TRUE)
        }
    }

    invisible(NULL)
}

#' Scale Values to \eqn{[0, 1]} Range
#'
#' Linearly scales values to the range \eqn{[0, 1]}
#'
#' @param x Numeric vector to scale.
#' @param low Minimum value for scaling. Defaults to min(x).
#' @param high Maximum value for scaling. Defaults to max(x).
#'
#' @return Numeric vector scaled to \eqn{[0, 1]} range.
#'
#' @examples
#' x <- c(1, 5, 10, 15, 20)
#' scale01(x)  # Returns: 0.00 0.21 0.47 0.74 1.00
#'
#' # Custom range
#' scale01(x, low = 0, high = 25)  # Different scaling
#'
#' @keywords internal
#' @noRd
scale01 <- function(x, low = min(x), high = max(x)) {

    if (!is.numeric(x)) {
        stop("x must be numeric")
    }

    if (low >= high) {
        stop("low must be less than high")
    }

    x <- (x - low) / (high - low)
    x
}

#' Plot Cluster Profile Lines
#'
#' Creates a line plot showing profiles of all members in a cluster
#'
#' @param X Matrix or data frame where rows are observations and columns are variables.
#' @param cltr Vector of cluster assignments.
#' @param id Cluster ID to plot.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#'
#' @return Invisibly returns the median profile.
#'
#' @details
#' This function creates a profile plot showing all members of a specified cluster
#' as gray lines, with the median profile highlighted in red. Missing values are
#' replaced with zeros.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), ncol = 10)
#' cltr <- sample(1:3, 10, replace = TRUE)
#' plot2D.cltr.profiles(X, cltr, id = 2, xlab = "Variables", ylab = "Values")
#' }
#'
#' @export
plot2D.cltr.profiles <- function(X,
                               cltr,
                               id,
                               xlab = "",
                               ylab = "") {

    if (length(cltr) != nrow(X)) {
        stop("Length of cltr must match number of rows in X")
    }

    # Identifying cluster profiles
    idx <- cltr == id

    if (sum(idx) == 0) {
        stop("No observations found for cluster id: ", id)
    }

    cl.profs <- X[idx, , drop = FALSE]

    # Replace NA's with 0s
    cl.profs <- apply(cl.profs, 2, function(x) {
        x[is.na(x)] <- 0
        x
    })

    # Ensure cl.profs is still a matrix
    if (!is.matrix(cl.profs)) {
        cl.profs <- matrix(cl.profs, nrow = 1)
    }

    # Line plots
    n <- ncol(cl.profs)
    ylim <- range(cl.profs, na.rm = TRUE)

    plot(1, 1, type = "n", xlim = c(1, n), ylim = ylim, las = 1,
         xlab = xlab, ylab = ylab)

    for (i in seq(nrow(cl.profs))) {
        lines(seq(n), cl.profs[i, ], col = "gray")
    }

    median.cl.profs <- apply(cl.profs, 2, median, na.rm = TRUE)
    lines(seq(n), median.cl.profs, col = "red", lwd = 3)

    invisible(median.cl.profs)
}

#' Plot Node Level Proportions on Graph
#'
#' Visualizes proportions of factor levels within graph nodes
#'
#' @param adj.mat Adjacency matrix of a graph.
#' @param props Matrix of proportions where rows are nodes and columns are factor levels.
#' @param color.tbl Vector of colors for factor levels. If NULL, colors are assigned automatically.
#' @param layout Character string specifying layout algorithm or a matrix of coordinates.
#'   Options include "auto", "circle", "sphere", "dh", "fr", "kk", "lgl", "mds",
#'   "sugiyama", "star", "tree", "grid", "random". Default is "auto".
#' @param epsilon Scaling factor for the radius of proportion circles within nodes.
#' @param legend.pos Position of the legend (e.g., "topleft", "topright").
#'
#' @return Invisibly returns the graph layout.
#'
#' @details
#' This function creates a network visualization where each node displays a pie chart
#' showing the proportions of different factor levels. The network structure is
#' determined by the adjacency matrix.
#'
#' @examples
#' \dontrun{
#' # Create example adjacency matrix
#' adj.mat <- matrix(0, 5, 5)
#' adj.mat\code{[1,2]} <- adj.mat\code{[2,1]} <- 1
#' adj.mat\code{[2,3]} <- adj.mat\code{[3,2]} <- 1
#'
#' # Create proportions matrix
#' props <- matrix(runif(15), ncol = 3)
#' props <- props / rowSums(props)
#' colnames(props) <- c("Type A", "Type B", "Type C")
#'
#' plot2D.node.level.props(adj.mat, props)
#' }
#'
#' @importFrom grDevices adjustcolor
#' @importFrom igraph graph_from_adjacency_matrix layout_nicely layout_in_circle
#' @importFrom igraph layout_on_sphere layout_with_dh layout_with_fr layout_with_kk
#' @importFrom igraph layout_with_lgl layout_with_mds layout_with_sugiyama
#' @importFrom igraph layout_as_star layout_as_tree layout_on_grid layout_randomly
#' @importFrom igraph plot.igraph
#' @export
plot2D.node.level.props <- function(adj.mat,
                                 props,
                                 color.tbl = NULL,
                                 layout = "auto",
                                 epsilon = 0.1,
                                 legend.pos = "topleft") {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required. Please install it.")
    }

    if (!is.matrix(adj.mat) || nrow(adj.mat) != ncol(adj.mat)) {
        stop("adj.mat must be a square matrix")
    }

    if (!is.matrix(props)) {
        stop("props must be a matrix")
    }

    if (nrow(adj.mat) != nrow(props)) {
        stop("Number of rows in props must match dimensions of adj.mat")
    }

    # Column names are level labels
    level.labels <- colnames(props)
    n.levels <- ncol(props)

    # Define colors for each level
    if (is.null(color.tbl)) {
        if (n.levels < 8) {
            color.tbl <- 2:(n.levels + 1)
        } else if (n.levels == 8) {
            color.tbl <- c(2:(n.levels + 1), "gray")
        } else {
            if (requireNamespace("mclust", quietly = TRUE)) {
                color.tbl <- mclust::mclust.options("classPlotColors")[seq(n.levels)]
            } else {
                color.tbl <- grDevices::rainbow(n.levels)
            }
        }
    }

    # Add transparency to colors
    color.tbl <- grDevices::adjustcolor(color.tbl, alpha.f = 0.5)

    # Create igraph object
    graph <- igraph::graph_from_adjacency_matrix(adj.mat, mode = "undirected",
                                                 weighted = NULL, diag = FALSE)

    # Calculate layout
    if (is.character(layout)) {
        graph_layout <- switch(layout,
            "auto" = igraph::layout_nicely(graph),
            "circle" = igraph::layout_in_circle(graph),
            "sphere" = igraph::layout_on_sphere(graph),
            "dh" = igraph::layout_with_dh(graph),
            "fr" = igraph::layout_with_fr(graph),
            "kk" = igraph::layout_with_kk(graph),
            "lgl" = igraph::layout_with_lgl(graph),
            "mds" = igraph::layout_with_mds(graph),
            "sugiyama" = igraph::layout_with_sugiyama(graph)$layout,
            "star" = igraph::layout_as_star(graph),
            "tree" = igraph::layout_as_tree(graph),
            "grid" = igraph::layout_on_grid(graph),
            "random" = igraph::layout_randomly(graph),
            igraph::layout_nicely(graph)  # default
        )
    } else if (is.matrix(layout)) {
        graph_layout <- layout
    } else {
        stop("layout must be a character string or a matrix")
    }

    # Scale layout to reasonable range
    graph_layout[,1] <- (graph_layout[,1] - min(graph_layout[,1])) /
                       (max(graph_layout[,1]) - min(graph_layout[,1])) * 4 - 2
    graph_layout[,2] <- (graph_layout[,2] - min(graph_layout[,2])) /
                       (max(graph_layout[,2]) - min(graph_layout[,2])) * 4 - 2

    # Plot the base graph
    plot(graph,
         layout = graph_layout,
         vertex.color = "white",
         vertex.size = 30,
         vertex.label = NA,  # We'll add labels manually
         edge.color = "gray",
         edge.width = 2,
         margin = 0.1)

    # Add colored pie segments for each node
    for (i in 1:nrow(props)) {
        xpos <- graph_layout[i, 1]
        ypos <- graph_layout[i, 2]

        start.angle <- 0
        for (j in 1:ncol(props)) {
            end.angle <- start.angle + props[i, j] * 2 * pi
            x.poly <- c(xpos, xpos + epsilon * cos(seq(start.angle, end.angle, length.out = 100)))
            y.poly <- c(ypos, ypos + epsilon * sin(seq(start.angle, end.angle, length.out = 100)))
            polygon(x.poly, y.poly, col = color.tbl[j], border = NA)
            start.angle <- end.angle
        }

        # Add node labels on top
        text(xpos, ypos, labels = if(!is.null(rownames(props))) rownames(props)[i] else i)
    }

    # Add legend
    legend(legend.pos, legend = level.labels, fill = color.tbl,
           inset = 0.05, cex = 0.8)

    invisible(graph_layout)
}

#' Draw Dashed Line in 3D Space
#'
#' Creates a dashed line between two points in 3D space
#'
#' @param x1,y1,z1 Coordinates of the starting point.
#' @param x2,y2,z2 Coordinates of the ending point.
#' @param n Deprecated parameter kept for compatibility.
#' @param dashLength Length of each dash.
#' @param gapLength Length of gaps between dashes.
#' @param col Color of the line.
#' @param lwd Line width.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function creates a dashed line by drawing multiple short line segments
#' with gaps between them.
#'
#' @examples
#' \dontrun{
#' # Draw a dashed line from origin to (1, 1, 1)
#' plot3D.plain(matrix(c(0,0,0,1,1,1), ncol = 3, byrow = TRUE))
#' draw.dashed.line3d(0, 0, 0, 1, 1, 1, col = "red", lwd = 3)
#' }
#'
#' @export
draw.dashed.line3d <- function(x1, y1, z1, x2, y2, z2,
                               n = 10,
                               dashLength = 0.1,
                               gapLength = 0.1,
                               col = "black",
                               lwd = 2) {

    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("Package 'rgl' is required. Please install it.")
    }

    # Calculate total length of the line segment
    totalLength <- sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)

    if (totalLength == 0) {
        warning("Start and end points are identical")
        return(invisible(NULL))
    }

    # Create sequence of distances along the line
    seqLengths <- seq(0, totalLength, by = dashLength + gapLength)

    for (i in seq_len(length(seqLengths) - 1)) {
        # Start and end points of the dashed segment
        startRatio <- seqLengths[i] / totalLength
        endRatio <- min(seqLengths[i] + dashLength, totalLength) / totalLength

        # Calculate coordinates
        startX <- x1 + (x2 - x1) * startRatio
        startY <- y1 + (y2 - y1) * startRatio
        startZ <- z1 + (z2 - z1) * startRatio
        endX <- x1 + (x2 - x1) * endRatio
        endY <- y1 + (y2 - y1) * endRatio
        endZ <- z1 + (z2 - z1) * endRatio

        # Draw the segment
        rgl::segments3d(c(startX, endX), c(startY, endY), c(startZ, endZ),
                       col = col, lwd = lwd)
    }

    invisible(NULL)
}

#' Create Circle Plot of Matrix Data
#'
#' Generates a circle plot where matrix columns are mapped to points on a unit circle
#'
#' @param X A data matrix where rows are observations and columns are variables.
#' @param n.circle.pts Number of points used to draw the circle.
#' @param delta Margin offset for accommodating labels.
#' @param pch Plotting character for data points.
#' @param cex Character expansion factor for data points.
#' @param col Color of data points.
#' @param comp.pch Plotting character for component points on circle.
#' @param comp.col Color of component points.
#' @param comp.cex Character expansion factor for component points.
#' @param lab.cex Character expansion factor for labels.
#' @param edge.col Color of edges from center to circle points.
#' @param adj.df Matrix of text adjustment values (n.cols x 2) for label positioning.
#'
#' @return Invisibly returns the component points matrix.
#'
#' @details
#' This function creates a visualization where the first variable is placed at the
#' center of a unit circle, and remaining variables are uniformly distributed around
#' the circle. Data points are then projected into this circular space.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(50), ncol = 5)
#' colnames(X) <- paste0("Var", 1:5)
#' circle.plot(X)
#' }
#'
#' @export
circle.plot <- function(X,
                        n.circle.pts = 100,
                        delta = 0.11,
                        pch = ".",
                        cex = 2.5,
                        col = "blue",
                        comp.pch = 19,
                        comp.col = "red",
                        comp.cex = 1.2,
                        lab.cex = 1.2,
                        edge.col = "gray",
                        adj.df = NULL) {

    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame")
    }

    p <- ncol(X)

    if (p < 2) {
        stop("X must have at least 2 columns")
    }

    # Generate circle points
    tt <- seq(0, 2 * pi, length = n.circle.pts)
    circle.df <- matrix(nrow = n.circle.pts, ncol = 2)
    for (i in seq(n.circle.pts)) {
        circle.df[i, ] <- c(cos(tt[i]), sin(tt[i]))
    }

    # Generate component points
    comp.pts <- matrix(0, nrow = p, ncol = 2)
    rownames(comp.pts) <- colnames(X)
    theta <- 2 * pi / (p - 1)
    for (i in 1:(p - 1)) {
       comp.pts[i + 1, ] <- c(cos(i * theta), sin(i * theta))
    }

    # Map X to circle
    X.circle <- as.matrix(X) %*% comp.pts

    # Create plot
    plot(circle.df, col = "gray", type = "l", axes = FALSE, xlab = "", ylab = "",
         xlim = c(-1 - delta, 1 + delta), ylim = c(-1 - delta, 1 + delta))
    box()

    # Draw edges
    for (i in 1:(p - 1)) {
        segments(x0 = 0, y0 = 0, x1 = comp.pts[i + 1, 1],
                y1 = comp.pts[i + 1, 2], col = edge.col)
    }

    # Add points
    points(X.circle, pch = pch, cex = cex, col = col)
    points(comp.pts, pch = comp.pch, col = comp.col, cex = comp.cex)

    # Add labels
    if (!is.null(adj.df)) {
        for (i in seq(p)) {
            text(comp.pts[i, 1], comp.pts[i, 2], labels = colnames(X)[i],
                 adj = adj.df[i, ], cex = lab.cex)
        }
    } else {
        for (i in seq(p)) {
            text(comp.pts[i, 1], comp.pts[i, 2], labels = colnames(X)[i],
                 cex = lab.cex)
        }
    }

    invisible(comp.pts)
}

#' Draw 3D Axes
#'
#' Creates and plots 3D axes with labels and arrows
#'
#' @param delta Extension length beyond unit length for axes.
#' @param axes.color Color of the axes.
#' @param axes.lab.cex Character expansion factor for axis labels.
#' @param edge.lwd Line width of axes.
#' @param half.axes Logical. If TRUE, only positive half-axes are drawn.
#' @param x.lab Label for x-axis.
#' @param y.lab Label for y-axis.
#' @param z.lab Label for z-axis.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function draws a set of 3D axes with optional arrows and labels.
#' The axes can be drawn as full axes or half-axes (positive direction only).
#'
#' @examples
#' \dontrun{
#' plot3D.plain(matrix(rnorm(30), ncol = 3))
#' draw.axes(delta = 0.5, axes.color = "black", half.axes = TRUE)
#' }
#'
#' @export
draw.axes <- function(delta = 0.5,
                      axes.color = "black",
                      axes.lab.cex = 2,
                      edge.lwd = 10,
                      half.axes = TRUE,
                      x.lab = "x_1",
                      y.lab = "x_2",
                      z.lab = "x_3") {

    axis.lim <- 1 + delta
    arrowLength <- 0.1

    if (half.axes) {
        # X-axis
        x.axis <- rbind(c(-delta, 0, 0), c(axis.lim, 0, 0))
        rgl::segments3d(x.axis, col = axes.color)
        rgl::text3d(x = axis.lim, y = 0, z = 0, texts = x.lab,
                   adj = c(1.5, 1.5), col = axes.color, cex = axes.lab.cex)

        # Y-axis
        y.axis <- rbind(c(0, -delta, 0), c(0, axis.lim, 0))
        rgl::segments3d(y.axis, col = axes.color)
        rgl::text3d(x = 0, y = axis.lim, z = 0, texts = y.lab,
                   adj = c(1.5, 1.5), col = axes.color, cex = axes.lab.cex)

        # Z-axis
        z.axis <- rbind(c(0, 0, -delta), c(0, 0, axis.lim))
        rgl::segments3d(z.axis, col = axes.color)
        rgl::text3d(x = 0, y = 0, z = axis.lim, texts = z.lab,
                   adj = c(1.5, 1.5), col = axes.color, cex = axes.lab.cex)
    } else {
        # Full axes
        x.axis <- rbind(c(-axis.lim, 0, 0), c(axis.lim, 0, 0))
        rgl::segments3d(x.axis, col = axes.color)
        rgl::text3d(x = axis.lim, y = 0, z = 0, texts = x.lab,
                   adj = c(1.5, 1.5), col = axes.color, cex = axes.lab.cex)

        y.axis <- rbind(c(0, -axis.lim, 0), c(0, axis.lim, 0))
        rgl::segments3d(y.axis, col = axes.color)
        rgl::text3d(x = 0, y = axis.lim, z = 0, texts = y.lab,
                   adj = c(1.5, 1.5), col = axes.color, cex = axes.lab.cex)

        z.axis <- rbind(c(0, 0, -axis.lim), c(0, 0, axis.lim))
        rgl::segments3d(z.axis, col = axes.color)
        rgl::text3d(x = 0, y = 0, z = axis.lim, texts = z.lab,
                   adj = c(1.5, 1.5), col = axes.color, cex = axes.lab.cex)
    }

    # Add arrows
    rgl::arrow3d(p0 = c(axis.lim - arrowLength, 0, 0),
                p1 = c(axis.lim, 0, 0),
                col = axes.color, type = "rotation", s = 1)
    rgl::arrow3d(p0 = c(0, axis.lim - arrowLength, 0),
                p1 = c(0, axis.lim, 0),
                col = axes.color, type = "rotation", s = 1)
    rgl::arrow3d(p0 = c(0, 0, axis.lim - arrowLength),
                p1 = c(0, 0, axis.lim),
                col = axes.color, type = "rotation", s = 1)

    invisible(NULL)
}

#' Draw 3D Line Segment
#'
#' Plots a line segment in 3D space from origin in the direction of a vector
#'
#' @param x Numeric vector of length 3 specifying direction.
#' @param length Length of the line segment.
#' @param col Color of the line segment.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function draws a line segment from the origin (0,0,0) in the direction
#' of vector x, with the specified length.
#'
#' @examples
#' \dontrun{
#' plot3D.plain(matrix(0, ncol = 3))
#' draw.3d.line(c(1, 1, 1), length = 2, col = "red")
#' draw.3d.line(c(1, 0, 0), length = 1.5, col = "blue")
#' }
#'
#' @export
draw.3d.line <- function(x, length = 2, col = "gray") {

    if (!is.numeric(x) || length(x) != 3) {
        stop("x must be a numeric vector of length 3")
    }

    if (!is.numeric(length) || length <= 0) {
        stop("length must be a positive number")
    }

    # Check for zero vector
    if (all(x == 0)) {
        stop("x cannot be a zero vector")
    }

    # Calculate L2 norm
    x.length <- sqrt(sum(x^2))

    # Rescale to desired length
    x <- length / x.length * x

    # Draw line segment
    rgl::segments3d(rbind(c(0, 0, 0), x), col = col)

    invisible(NULL)
}

#' Plot 3d Disk Embedding Object
#'
#' Visualizes a disk embedding object in 3D space
#'
#' @param ebdg.obj A diskEmbdg object containing embedding information.
#' @param col Color of spheres representing data points.
#' @param radius Radius of data point spheres.
#' @param vertex.radius Radius of vertex spheres.
#' @param vertex.col Color of vertex spheres.
#' @param vertex.cex Character expansion factor for vertex labels.
#' @param edge.col Color of edges.
#' @param adj.df Matrix of text adjustment values for vertex labels.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function visualizes a disk embedding where data is projected onto a
#' sphere with axes represented as points on the sphere and one at the center.
#'
#' @examples
#' \dontrun{
#' # Assuming ebdg.obj is a properly formatted disk embedding object
#' plot3D.diskEmbdg(ebdg.obj, col = "blue", vertex.col = "red")
#' }
#'
#' @export
plot3D.diskEmbdg <- function(ebdg.obj,
                           col = "blue",
                           radius = 0.005,
                           vertex.radius = 0.02,
                           vertex.col = "red",
                           vertex.cex = 1.2,
                           edge.col = "gray",
                           adj.df = NULL) {

    # Extract components
    axis.pos <- ebdg.obj$axis.pos
    X.ebdg <- ebdg.obj$X.ebdg
    X <- ebdg.obj$X
    n <- ebdg.obj$n

    # Create sphere plot
    plot3D.plain(X.ebdg)
    rgl::spheres3d(X.ebdg, col = col, radius = radius)

    # Draw line segments from center to sphere points
    for (i in seq(n)) {
        rgl::segments3d(rbind(numeric(3), axis.pos[i, ]), col = edge.col)
    }

    # Draw vertices
    rgl::spheres3d(axis.pos, col = vertex.col, radius = vertex.radius)

    # Draw labels
    if (!is.null(adj.df)) {
        for (i in seq(n)) {
            rgl::text3d(axis.pos[i, ], texts = colnames(X)[i],
                       cex = vertex.cex, adj = adj.df[i, ])
        }
    } else {
        for (i in seq(n)) {
            rgl::text3d(axis.pos[i, ], texts = colnames(X)[i], cex = vertex.cex)
        }
    }

    invisible(NULL)
}

#' Add Vertical Error Bar to Plot
#'
#' Adds a vertical error bar to an existing 2D plot
#'
#' @param x X-coordinate for the error bar.
#' @param ymin Lower limit of error bar.
#' @param ymax Upper limit of error bar.
#' @param dx Horizontal offset for end caps.
#' @param lwd Line width.
#' @param col Color of error bar.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' plot(1:10, rnorm(10), ylim = c(-3, 3))
#' vert.error.bar(5, -1, 1, col = "red")
#'
#' @export
vert.error.bar <- function(x, ymin, ymax, dx = 0.025, lwd = 1, col = "red") {
    segments(x, ymin, x, ymax, col = col, lwd = lwd)
    segments(x - dx, ymin, x + dx, ymin, col = col, lwd = lwd)
    segments(x - dx, ymax, x + dx, ymax, col = col, lwd = lwd)
    invisible(NULL)
}

#' Add Horizontal Error Bar to Plot
#'
#' Adds a horizontal error bar to an existing 2D plot
#'
#' @param xmin Left limit of error bar.
#' @param xmax Right limit of error bar.
#' @param y Y-coordinate for the error bar.
#' @param dyy Vertical offset for end caps.
#' @param lwd Line width.
#' @param col Color of error bar.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' plot(rnorm(10), 1:10, xlim = c(-3, 3))
#' hor.error.bar(-1, 1, 5, col = "blue")
#'
#' @export
hor.error.bar <- function(xmin, xmax, y, dyy = 0.025, lwd = 1, col = "black") {
    segments(xmin, y, xmax, y, col = col, lwd = lwd)
    segments(xmin, y - dyy, xmin, y + dyy, col = col, lwd = lwd)
    segments(xmax, y - dyy, xmax, y + dyy, col = col, lwd = lwd)
    invisible(NULL)
}

#' Panel Function for Absolute Correlations
#'
#' Displays absolute correlations in scatter plot matrix panels
#'
#' @param x Numeric vector for first variable.
#' @param y Numeric vector for second variable.
#' @param digits Number of digits to display.
#' @param prefix Text prefix for correlation value.
#' @param cex.cor Base character expansion factor.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This function is designed to be used with the pairs() function to show
#' absolute correlations in the upper or lower panels. Text size is proportional
#' to the correlation magnitude.
#'
#' @examples
#' pairs(iris[1:4], upper.panel = panel.acor)
#'
#' @keywords internal
#' @noRd
panel.acor <- function(x, y, digits = 2, prefix = "", cex.cor = 1, ...) {

    # Handle missing values
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]

    if (length(x) < 3) {
        return(invisible(NULL))
    }

    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))

    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)

    if (missing(cex.cor)) {
        cex.cor <- 0.8 / strwidth(txt)
    }

    text(0.5, 0.5, txt, cex = cex.cor * r)

    invisible(NULL)
}

#' Panel Function with Robust Line
#'
#' Adds scatter plot with Tukey's robust line to pairs() panels
#'
#' @param x Numeric vector for first variable.
#' @param y Numeric vector for second variable.
#' @param col Color for points.
#' @param bg Background color for points.
#' @param pch Plotting character.
#' @param cex Character expansion factor for points.
#' @param col.rlm Color for robust line.
#' @param ... Additional parameters passed to abline().
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' This panel function creates a scatter plot with Tukey's resistant line
#' fitted to the data. Useful for pairs() plots to show robust relationships.
#'
#' @examples
#' pairs(iris[1:4], lower.panel = panel.rlm)
#'
#' @keywords internal
#' @noRd
panel.rlm <- function(x, y, col = par("col"), bg = NA, pch = par("pch"),
                     cex = 1, col.rlm = 2, ...) {

    points(x, y, pch = pch, col = col, bg = bg, cex = cex)

    ok <- is.finite(x) & is.finite(y)
    if (any(ok) && length(x[ok]) > 2) {
        abline(line(x[ok], y[ok]), col = col.rlm, ...)
    }

    invisible(NULL)
}


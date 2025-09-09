#' Histogram with customized default settings
#'
#' Creates a histogram with customized default settings for n.breaks, color, and appearance.
#'
#' @param x A numeric vector of values for which the histogram is plotted.
#' @param main The title of the plot. Default is an empty string.
#' @param n.breaks The number of breaks (bins) in the histogram. Default is 100.
#' @param col The color of the histogram bars. Default is "red".
#' @param ... Additional parameters passed to \code{\link[graphics]{hist}}.
#'
#' @return Invisibly returns a list of class "histogram" as returned by \code{\link[graphics]{hist}}.
#'
#' @export
#' @importFrom graphics hist
#'
#' @examples
#' # Basic usage
#' data <- rnorm(1000)
#' hist1(data)
#'
#' # With custom settings
#' hist1(data, main = "Normal Distribution", n.breaks = 50, col = "blue")
#'
#' # Using additional parameters
#' hist1(data, main = "Custom Histogram", n.breaks = 30, col = "green",
#'       xlab = "Values", ylab = "Frequency")
#'
hist1 <- function(x, main = "", n.breaks = 100, col = "red", ...) {
    # Input validation
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector")
    }
    if (length(x) == 0) {
        stop("'x' must contain at least one value")
    }
    if (!is.numeric(n.breaks) || n.breaks <= 0) {
        stop("'n.breaks' must be a positive number")
    }

    # Create histogram
    invisible(hist(x, breaks = n.breaks, las = 1, main = main, col = col, ...))
}


#' Plot two overlapping histograms
#'
#' Creates two overlapping histograms with customizable colors and transparency
#' for comparing two distributions.
#'
#' @param x1 A numeric vector for the first histogram.
#' @param x2 A numeric vector for the second histogram.
#' @param x1.lab Label for x1 in the legend.
#' @param x2.lab Label for x2 in the legend.
#' @param xlab Label for the x-axis. Default is an empty string.
#' @param n.x1.breaks The number of breaks for the first histogram. Default is 30.
#' @param n.x2.breaks The number of breaks for the second histogram. Default is 30.
#' @param main The title of the plot. Default is an empty string.
#' @param col.x1 Color of x1 histogram bars. Default is light blue with transparency.
#' @param col.x2 Color of x2 histogram bars. Default is light pink with transparency.
#' @param legend.pos Position of the legend. Default is "topright".
#'              Must be one of "bottomright", "bottom", "bottomleft", "left",
#'              "topleft", "top", "topright", "right", "center" or "none".
#'
#' @return Invisibly returns a list containing:
#'   \item{xlim}{The x-axis limits used in the plot}
#'   \item{ylim}{The y-axis limits used in the plot}
#'
#' @export
#' @importFrom graphics hist legend
#' @importFrom grDevices rgb
#'
#' @examples
#' # Basic usage
#' x1 <- rnorm(1000, mean = 0, sd = 1)
#' x2 <- rnorm(1000, mean = 2, sd = 1.5)
#' hist2(x1, x2, x1.lab = "Group A", x2.lab = "Group B")
#'
#' # With custom settings
#' hist2(x1, x2, x1.lab = "Control", x2.lab = "Treatment",
#'       xlab = "Values", main = "Comparison of Distributions",
#'       n.x1.breaks = 50, n.x2.breaks = 50)
#'
#' # Custom colors and no legend
#' hist2(x1, x2, x1.lab = "A", x2.lab = "B",
#'       col.x1 = "lightgreen", col.x2 = "lightcoral",
#'       legend.pos = "none")
#'
hist2 <- function(x1,
                  x2,
                  x1.lab,
                  x2.lab,
                  xlab = "",
                  n.x1.breaks = 30,
                  n.x2.breaks = 30,
                  main = "",
                  col.x1 = rgb(173, 216, 230, maxColorValue = 255, alpha = 80),
                  col.x2 = rgb(255, 192, 203, maxColorValue = 255, alpha = 80),
                  legend.pos = "topright") {

    # Input validation
    if (!is.numeric(x1) || !is.numeric(x2)) {
        stop("'x1' and 'x2' must be numeric vectors")
    }
    if (length(x1) == 0 || length(x2) == 0) {
        stop("'x1' and 'x2' must contain at least one value")
    }
    if (missing(x1.lab) || missing(x2.lab)) {
        stop("'x1.lab' and 'x2.lab' must be provided")
    }
    if (!is.numeric(n.x1.breaks) || !is.numeric(n.x2.breaks) ||
        n.x1.breaks <= 0 || n.x2.breaks <= 0) {
        stop("'n.x1.breaks' and 'n.x2.breaks' must be positive numbers")
    }

    valid.positions <- c("bottomright", "bottom", "bottomleft", "left",
                        "topleft", "top", "topright", "right", "center", "none")
    if (!legend.pos %in% valid.positions) {
        stop("'legend.pos' must be one of: ", paste(valid.positions, collapse = ", "))
    }

    # Calculate densities to determine y-axis limits
    hist.x1 <- hist(x1, breaks = n.x1.breaks, plot = FALSE)
    hist.x2 <- hist(x2, breaks = n.x2.breaks, plot = FALSE)

    max.x1.den <- max(hist.x1$density)
    max.x2.den <- max(hist.x2$density)

    # Set plot limits
    xlim <- range(c(x1, x2))
    ylim <- c(0, max(c(max.x1.den, max.x2.den)))

    # Create first histogram
    hist(x1,
         breaks = n.x1.breaks,
         probability = TRUE,
         col = col.x1,
         las = 1,
         xlab = xlab,
         xlim = xlim,
         ylim = ylim,
         main = main)

    # Add second histogram
    hist(x2, breaks = n.x2.breaks, probability = TRUE, col = col.x2, add = TRUE)

    # Add legend if requested
    if (legend.pos != "none") {
        legend(legend.pos,
               legend = c(x1.lab, x2.lab),
               fill = c(col.x1, col.x2),
               bty = "n",
               inset = 0.01)
    }

    invisible(list(xlim = xlim, ylim = ylim))
}


#' Plot three overlapping histograms
#'
#' Creates three overlapping histograms with customizable colors and transparency
#' for comparing three distributions.
#'
#' @param x1 A numeric vector for the first histogram.
#' @param x2 A numeric vector for the second histogram.
#' @param x3 A numeric vector for the third histogram.
#' @param x1.lab Label for x1 in the legend.
#' @param x2.lab Label for x2 in the legend.
#' @param x3.lab Label for x3 in the legend.
#' @param xlab Label for the x-axis. Default is an empty string.
#' @param n.x1.breaks The number of breaks for the first histogram. Default is 30.
#' @param n.x2.breaks The number of breaks for the second histogram. Default is 30.
#' @param n.x3.breaks The number of breaks for the third histogram. Default is 30.
#' @param main The title of the plot. Default is an empty string.
#' @param ylim.max The upper y-axis limit. If NA (default), automatically determined.
#' @param col.x1 Color of x1 histogram bars. Default is semi-transparent red.
#' @param col.x2 Color of x2 histogram bars. Default is semi-transparent gray.
#' @param col.x3 Color of x3 histogram bars. Default is semi-transparent blue.
#' @param legend.pos Position of the legend. Default is "topright".
#'              Must be one of "bottomright", "bottom", "bottomleft", "left",
#'              "topleft", "top", "topright", "right", "center" or "none".
#'
#' @return Invisibly returns a list containing:
#'   \item{xlim}{The x-axis limits used in the plot}
#'   \item{ylim}{The y-axis limits used in the plot}
#'
#' @export
#' @importFrom graphics hist legend
#' @importFrom grDevices rgb
#'
#' @examples
#' # Basic usage
#' x1 <- rnorm(1000, mean = 0, sd = 1)
#' x2 <- rnorm(1000, mean = 2, sd = 1.2)
#' x3 <- rnorm(1000, mean = 1, sd = 0.8)
#' hist3(x1, x2, x3,
#'       x1.lab = "Group A", x2.lab = "Group B", x3.lab = "Group C")
#'
#' # With custom settings
#' hist3(x1, x2, x3,
#'       x1.lab = "Control", x2.lab = "Treatment 1", x3.lab = "Treatment 2",
#'       xlab = "Measurements", main = "Three-way Comparison",
#'       n.x1.breaks = 50, n.x2.breaks = 50, n.x3.breaks = 50)
#'
#' # Fixed y-axis limit and custom colors
#' hist3(x1, x2, x3,
#'       x1.lab = "A", x2.lab = "B", x3.lab = "C",
#'       ylim.max = 0.5,
#'       col.x1 = "lightgreen", col.x2 = "lightblue", col.x3 = "lightcoral")
#'
hist3 <- function(x1, x2, x3,
                  x1.lab, x2.lab, x3.lab,
                  xlab = "",
                  n.x1.breaks = 30,
                  n.x2.breaks = 30,
                  n.x3.breaks = 30,
                  main = "",
                  ylim.max = NA,
                  col.x1 = rgb(255, 0, 0, maxColorValue = 255, alpha = 90),
                  col.x2 = rgb(190, 190, 190, maxColorValue = 255, alpha = 90),
                  col.x3 = rgb(0, 0, 255, maxColorValue = 255, alpha = 90),
                  legend.pos = "topright") {

    # Input validation
    if (!is.numeric(x1) || !is.numeric(x2) || !is.numeric(x3)) {
        stop("'x1', 'x2', and 'x3' must be numeric vectors")
    }
    if (length(x1) == 0 || length(x2) == 0 || length(x3) == 0) {
        stop("'x1', 'x2', and 'x3' must contain at least one value")
    }
    if (missing(x1.lab) || missing(x2.lab) || missing(x3.lab)) {
        stop("'x1.lab', 'x2.lab', and 'x3.lab' must be provided")
    }
    if (!is.numeric(n.x1.breaks) || !is.numeric(n.x2.breaks) || !is.numeric(n.x3.breaks) ||
        n.x1.breaks <= 0 || n.x2.breaks <= 0 || n.x3.breaks <= 0) {
        stop("'n.x1.breaks', 'n.x2.breaks', and 'n.x3.breaks' must be positive numbers")
    }
    if (!is.na(ylim.max) && (!is.numeric(ylim.max) || ylim.max <= 0)) {
        stop("'ylim.max' must be NA or a positive number")
    }

    valid.positions <- c("bottomright", "bottom", "bottomleft", "left",
                        "topleft", "top", "topright", "right", "center", "none")
    if (!legend.pos %in% valid.positions) {
        stop("'legend.pos' must be one of: ", paste(valid.positions, collapse = ", "))
    }

    # Calculate densities to determine y-axis limits
    hist.x1 <- hist(x1, breaks = n.x1.breaks, plot = FALSE)
    hist.x2 <- hist(x2, breaks = n.x2.breaks, plot = FALSE)
    hist.x3 <- hist(x3, breaks = n.x3.breaks, plot = FALSE)

    max.x1.den <- max(hist.x1$density)
    max.x2.den <- max(hist.x2$density)
    max.x3.den <- max(hist.x3$density)

    # Set plot limits
    xlim <- range(c(x1, x2, x3))

    if (is.na(ylim.max)) {
        ylim.max <- max(c(max.x1.den, max.x2.den, max.x3.den))
    }
    ylim <- c(0, ylim.max)

    # Create first histogram
    hist(x1,
         breaks = n.x1.breaks,
         probability = TRUE,
         col = col.x1,
         xlab = xlab,
         xlim = xlim,
         las = 1,
         ylim = ylim,
         main = main)

    # Add second and third histograms
    hist(x2, breaks = n.x2.breaks, probability = TRUE, col = col.x2, add = TRUE)
    hist(x3, breaks = n.x3.breaks, probability = TRUE, col = col.x3, add = TRUE)

    # Add legend if requested
    if (legend.pos != "none") {
        legend(legend.pos,
               legend = c(x1.lab, x2.lab, x3.lab),
               fill = c(col.x1, col.x2, col.x3),
               bty = "n",
               inset = 0.01)
    }

    invisible(list(xlim = xlim, ylim = ylim))
}

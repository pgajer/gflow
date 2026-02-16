#' Apply Sigmoid Transformation
#'
#' Transforms input values using one of three sigmoid functions that map
#' to the range (-1, 1). These transformations are useful for normalizing
#' correlation-like quantities while preserving sign and relative magnitude.
#'
#' @param x Numeric vector of input values.
#' @param alpha Positive numeric scale parameter controlling the steepness
#'   of the transformation. Larger values produce steeper transitions near zero.
#' @param type Character string specifying the sigmoid type. One of:
#'   \describe{
#'     \item{\code{"tanh"}}{Hyperbolic tangent: \eqn{\tanh(\alpha x)}}
#'     \item{\code{"arctan"}}{Scaled arctangent: \eqn{(2/\pi) \arctan(\alpha x)}}
#'     \item{\code{"algebraic"}}{Algebraic sigmoid: \eqn{x / \sqrt{1/\alpha^2 + x^2}}}
#'   }
#'
#' @return Numeric vector of transformed values in (-1, 1).
#'
#' @details
#' All three sigmoid functions share key properties: they are odd functions
#' (antisymmetric about the origin), monotonically increasing, and bounded
#' by (-1, 1). The parameter \code{alpha} controls the rate of saturation,
#' with larger values causing faster approach to the asymptotic bounds.
#'
#' The hyperbolic tangent is the most common choice in practice. The arctangent
#' variant has slightly lighter tails. The algebraic form avoids transcendental
#' functions entirely, which can offer computational advantages in some contexts.
#'
#' @examples
#' x <- seq(-3, 3, length.out = 101)
#' y.tanh <- apply.sigmoid(x, alpha = 1, type = "tanh")
#' y.arctan <- apply.sigmoid(x, alpha = 1, type = "arctan")
#' y.algebraic <- apply.sigmoid(x, alpha = 1, type = "algebraic")
#'
#' @export
apply.sigmoid <- function(x, alpha = 1.0, type = c("tanh", "arctan", "algebraic")) {
    ## Validate inputs
    if (!is.numeric(x)) {
        stop("'x' must be numeric")
    }
    if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
        stop("'alpha' must be a single positive number")
    }

    type <- match.arg(type)

    ## Apply the appropriate transformation
    result <- switch(type,
        tanh = tanh(alpha * x),
        arctan = (2 / pi) * atan(alpha * x),
        algebraic = {
            inv.alpha.sq <- 1 / (alpha * alpha)
            x / sqrt(inv.alpha.sq + x * x)
        }
    )

    return(result)
}


#' Plot Comparison of Sigmoid Transformations
#'
#' Creates a grid of plots comparing the three sigmoid transformation types
#' across different values of the scale parameter alpha.
#'
#' @param alpha.values Numeric vector of alpha values to display. Each value
#'   generates one panel in the plot grid.
#' @param x.range Numeric vector of length 2 specifying the range of x values
#'   to plot. Default is \code{c(-3, 3)}.
#' @param n.points Integer number of points for smooth curve rendering.
#'   Default is 201.
#' @param cols Character vector of length 3 specifying colors for tanh,
#'   arctan, and algebraic curves respectively.
#' @param lwd Line width for curves. Default is 2.
#' @param add.legend Logical indicating whether to add a legend to each panel.
#'   Default is \code{TRUE}.
#' @param main.title Optional character string for an overall title.
#'
#' @return Invisible NULL. Called for its side effect of producing a plot.
#'
#' @examples
#' ## Default 3x3 comparison
#' plot.sigmoid.comparison()
#'
#' ## Custom alpha values
#' plot.sigmoid.comparison(alpha.values = c(0.25, 0.5, 1, 2, 4, 8))
#'
#' @export
plot.sigmoid.comparison <- function(alpha.values = c(0.5, 1, 2, 3, 5, 8, 12, 20, 50),
                                    x.range = c(-3, 3),
                                    n.points = 201,
                                    cols = c("#2E86AB", "#A23B72", "#F18F01"),
                                    lwd = 2,
                                    add.legend = TRUE,
                                    main.title = NULL) {
    ## Validate inputs
    if (!is.numeric(alpha.values) || length(alpha.values) < 1) {
        stop("'alpha.values' must be a numeric vector with at least one element")
    }
    if (any(alpha.values <= 0)) {
        stop("All alpha values must be positive")
    }

    n.panels <- length(alpha.values)

    ## Determine grid dimensions
    n.cols <- ceiling(sqrt(n.panels))
    n.rows <- ceiling(n.panels / n.cols)

    ## Set up the plotting grid
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    if (!is.null(main.title)) {
        par(mfrow = c(n.rows, n.cols),
            mar = c(4, 4, 3, 1),
            oma = c(0, 0, 3, 0))
    } else {
        par(mfrow = c(n.rows, n.cols),
            mar = c(4, 4, 3, 1))
    }

    ## Generate x values
    x <- seq(x.range[1], x.range[2], length.out = n.points)

    ## Sigmoid type labels for legend
    type.labels <- c("tanh", "arctan", "algebraic")

    ## Create each panel
    for (i in seq_along(alpha.values)) {
        alpha <- alpha.values[i]

        ## Compute all three sigmoid curves
        y.tanh <- apply.sigmoid(x, alpha = alpha, type = "tanh")
        y.arctan <- apply.sigmoid(x, alpha = alpha, type = "arctan")
        y.algebraic <- apply.sigmoid(x, alpha = alpha, type = "algebraic")

        ## Create the plot
        plot(x, y.tanh, type = "l", col = cols[1], lwd = lwd,
             xlim = x.range, ylim = c(-1.1, 1.1),
             xlab = "x", ylab = expression(sigma(x)),
             main = bquote(alpha == .(alpha)))

        ## Add other curves
        lines(x, y.arctan, col = cols[2], lwd = lwd, lty = 2)
        lines(x, y.algebraic, col = cols[3], lwd = lwd, lty = 3)

        ## Add reference lines
        abline(h = c(-1, 0, 1), col = "gray70", lty = 3, lwd = 0.5)
        abline(v = 0, col = "gray70", lty = 3, lwd = 0.5)

        ## Add legend to first panel only (or all if desired)
        if (add.legend && i == 1) {
            legend("bottomright", legend = type.labels,
                   col = cols, lty = c(1, 2, 3), lwd = lwd,
                   cex = 0.7, bg = "white")
        }
    }

    ## Add overall title if specified
    if (!is.null(main.title)) {
        mtext(main.title, outer = TRUE, cex = 1.2, line = 1)
    }

    invisible(NULL)
}


#' Plot Sigmoid Sensitivity to Alpha
#'
#' Creates a single plot showing how one sigmoid type varies with different
#' alpha values, useful for understanding the effect of the scale parameter.
#'
#' @param type Character string specifying the sigmoid type.
#' @param alpha.values Numeric vector of alpha values to compare.
#' @param x.range Numeric vector of length 2 specifying the x range.
#' @param n.points Integer number of points for curve rendering.
#' @param col.palette Function that generates colors, or a character vector.
#'
#' @return Invisible NULL.
#'
#' @examples
#' plot.sigmoid.alpha.sensitivity(type = "tanh",
#'                                alpha.values = c(0.5, 1, 2, 5, 10))
#'
#' @export
plot.sigmoid.alpha.sensitivity <- function(type = c("tanh", "arctan", "algebraic"),
                                           alpha.values = c(0.5, 1, 2, 5, 10),
                                           x.range = c(-3, 3),
                                           n.points = 201,
                                           col.palette = NULL) {
    type <- match.arg(type)
    n.alpha <- length(alpha.values)

    ## Generate colors
    if (is.null(col.palette)) {
        cols <- hcl.colors(n.alpha, palette = "viridis")
    } else if (is.function(col.palette)) {
        cols <- col.palette(n.alpha)
    } else {
        cols <- rep_len(col.palette, n.alpha)
    }

    ## Generate x values
    x <- seq(x.range[1], x.range[2], length.out = n.points)

    ## Set up the plot
    plot(NULL, xlim = x.range, ylim = c(-1.1, 1.1),
         xlab = "x", ylab = expression(sigma(x)),
         main = paste0("Sigmoid type: ", type))

    ## Add reference lines
    abline(h = c(-1, 0, 1), col = "gray80", lty = 3)
    abline(v = 0, col = "gray80", lty = 3)

    ## Add curves for each alpha
    for (i in seq_along(alpha.values)) {
        y <- apply.sigmoid(x, alpha = alpha.values[i], type = type)
        lines(x, y, col = cols[i], lwd = 2)
    }

    ## Add legend
    legend("bottomright",
           legend = paste0(expression(alpha), " = ", alpha.values),
           col = cols, lwd = 2, cex = 0.8, bg = "white")

    invisible(NULL)
}

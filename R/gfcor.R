#' Gradient Flow Correlation Analysis
#'
#' @description
#' Computes symmetric association measures between two fitted surfaces using
#' gradient flow complex structure. The analysis quantifies how the gradient
#' flow partitions of two functions relate to each other, providing measures
#' of association at global, basin, and vertex levels.
#'
#' @details
#' The gradient flow correlation framework addresses a fundamental question:
#' given two fitted surfaces y and z on a graph, in what regions do they
#' exhibit positive association (both increase together) versus negative
#' association (one increases while the other decreases)?
#'
#' The key insight is that gradient flow partitions the graph into basins of
#' attraction. Each basin consists of vertices whose gradient flow trajectories
#' converge to the same local extremum. By analyzing how the basin partitions
#' of y and z overlap and how the functions vary within those overlapping
#' regions, we obtain robust association measures that are insensitive to
#' spurious local extrema.
#'
#' The polarity coordinate p(v) in \eqn{[-1, 1]} measures where vertex v sits within
#' its accessible dynamic range. A vertex with p(v) near +1 sits close to a
#' local maximum of its gradient flow cell, while p(v) near -1 indicates
#' proximity to a local minimum. The vertex-level association a_pol(v) is the
#' product \eqn{p_y(v) \cdot p_z(v)}, which is positive when v occupies similar relative
#' positions in both landscapes.
#'
#' Two global summaries aggregate vertex-level association. The polarity
#' concordance A_pol is the mass-weighted average of a_pol(v), while the sign
#' concordance kappa_pol counts the mass-weighted proportion of vertices with
#' positive minus negative association signs.
#'
#' Basin association character chi provides a basin-level summary. For each
#' y-basin, chi measures the average z-polarity within that basin. Values near
#' +1 indicate that high-y regions coincide with high-z regions, while values
#' near -1 indicate opposite association.
#'
#' @section Polarity Scale:
#' Two polarity computation modes are available. The value-based mode (default)
#' uses raw function values to compute normalized height within cells. The
#' rank-based mode uses ranks within each cell, providing invariance under
#' monotone transformations of the functions. Rank-based polarity is recommended
#' for cross-cohort comparisons where absolute scales may differ.
#'
#' @section Basin Multiplicity:
#' In discrete Morse theory, gradient trajectories are not unique. A vertex may
#' belong to multiple basins simultaneously when it has multiple neighbors with
#' lower (or higher) function values. The implementation handles this multiplicity
#' through soft membership vectors that distribute unit mass across containing
#' basins. All statistics properly account for this multiplicity.
#'
#' @param y.hat Numeric vector of fitted values for the first surface. Length
#'   must equal the number of vertices in the graph.
#' @param z.hat Numeric vector of fitted values for the second surface. Must
#'   have the same length as y.hat.
#' @param y.basins Object of class \code{"basins_of_attraction"} computed from
#'   y.hat using \code{compute.basins.of.attraction}.
#' @param z.basins Object of class \code{"basins_of_attraction"} computed from
#'   z.hat using \code{compute.basins.of.attraction}.
#' @param vertex.mass Optional numeric vector of vertex weights for computing
#'   weighted averages. If NULL (default), uniform weights are used. Useful for
#'   incorporating sampling density or confidence weights.
#' @param polarity.scale Character string specifying the polarity computation
#'   mode. Either \code{"value"} (default) for normalized height based on
#'   function values, or \code{"rank"} for rank-based computation within cells.
#' @param epsilon Numeric threshold for detecting flat regions where polarity
#'   is undefined. Vertices with dynamic range smaller than epsilon are marked
#'   as invalid. Default is 1e-10.
#'
#' @return An object of class \code{"gfcor"} containing:
#'   \item{global}{List of global association measures:
#'     \itemize{
#'       \item \code{A_pol}: Polarity concordance in \eqn{[-1, 1]}
#'       \item \code{kappa_pol}: Sign concordance in \eqn{[-1, 1]}
#'       \item \code{n_positive}: Count of vertices with positive association
#'       \item \code{n_negative}: Count of vertices with negative association
#'       \item \code{n_zero}: Count of vertices with zero association
#'       \item \code{n_invalid}: Count of vertices with undefined polarity
#'     }}
#'   \item{vertex}{List of vertex-level results:
#'     \itemize{
#'       \item \code{a_pol}: Association score p_y * p_z in \eqn{[-1, 1]}
#'       \item \code{sign_pol}: Sign of association in \eqn{\{-1, 0, +1\}}
#'       \item \code{confidence}: Absolute value |a_pol| as confidence proxy
#'       \item \code{is_valid}: Logical indicating valid polarity
#'     }}
#'   \item{polarity_y}{Polarity structure for y (theta, polarity, range, is_valid)}
#'   \item{polarity_z}{Polarity structure for z}
#'   \item{basin_character}{Basin association characters:
#'     \itemize{
#'       \item \code{chi_y_max}: Character of each y-maximum basin
#'       \item \code{chi_y_min}: Character of each y-minimum basin
#'       \item \code{chi_z_max}: Character of each z-maximum basin
#'       \item \code{chi_z_min}: Character of each z-minimum basin
#'       \item \code{mass_*}: Corresponding basin masses
#'     }}
#'   \item{overlap}{Soft overlap matrices:
#'     \itemize{
#'       \item \code{O_pp}: y-max with z-max overlap
#'       \item \code{O_mm}: y-min with z-min overlap
#'       \item \code{O_pm}: y-max with z-min overlap
#'       \item \code{O_mp}: y-min with z-max overlap
#'       \item \code{total_mass}: Sum of vertex masses
#'     }}
#'   \item{membership}{List containing y and z membership structures}
#'
#' @examples
#' \dontrun{
#' ## Compute basins for two fitted surfaces
#' y.basins <- compute.basins.of.attraction(adj.list, weight.list, y.hat)
#' z.basins <- compute.basins.of.attraction(adj.list, weight.list, z.hat)
#'
#' ## Compute gradient flow correlation
#' gfc <- gfcor(y.hat, z.hat, y.basins, z.basins)
#'
#' ## Examine global association
#' print(gfc$global$A_pol)      ## Polarity concordance
#' print(gfc$global$kappa_pol)  ## Sign concordance
#'
#' ## Map vertex-level association
#' assoc.map <- gfc$vertex$a_pol
#' positive.vertices <- which(gfc$vertex$sign_pol > 0 & gfc$vertex$is_valid)
#'
#' ## Examine basin characters
#' print(gfc$basin_character$chi_y_max)  ## Which y-max basins are positively associated?
#' }
#'
#' @seealso
#' \code{\link{compute.basins.of.attraction}} for computing gradient flow basins,
#' \code{\link{gfcomon}} for directed co-monotonicity analysis,
#' \code{\link{gfassoc}} for comprehensive association analysis
#'
#' @export
gfcor <- function(y.hat,
                  z.hat,
                  y.basins,
                  z.basins,
                  vertex.mass = NULL,
                  polarity.scale = c("value", "rank"),
                  epsilon = 1e-10) {

    ## Input validation
    polarity.scale <- match.arg(polarity.scale)

    if (!is.numeric(y.hat) || !is.numeric(z.hat)) {
        stop("y.hat and z.hat must be numeric vectors")
    }

    if (length(y.hat) != length(z.hat)) {
        stop("y.hat and z.hat must have the same length")
    }

    if (!inherits(y.basins, "basins_of_attraction")) {
        stop("y.basins must be of class 'basins_of_attraction'")
    }

    if (!inherits(z.basins, "basins_of_attraction")) {
        stop("z.basins must be of class 'basins_of_attraction'")
    }

    if (y.basins$n_vertices != length(y.hat)) {
        stop("y.basins$n_vertices must match length of y.hat")
    }

    if (z.basins$n_vertices != length(z.hat)) {
        stop("z.basins$n_vertices must match length of z.hat")
    }

    if (!is.null(vertex.mass)) {
        if (!is.numeric(vertex.mass) || length(vertex.mass) != length(y.hat)) {
            stop("vertex.mass must be a numeric vector of the same length as y.hat")
        }
        if (any(vertex.mass < 0)) {
            stop("vertex.mass must be non-negative")
        }
    }

    if (!is.numeric(epsilon) || length(epsilon) != 1 || epsilon <= 0) {
        stop("epsilon must be a positive numeric scalar")
    }

    ## Build options list
    options <- list(
        polarity.scale = polarity.scale,
        epsilon = epsilon
    )

    ## Call C++ implementation
    result <- .Call(
        S_gfcor,
        as.numeric(y.hat),
        as.numeric(z.hat),
        y.basins$lmax_basins,
        y.basins$lmin_basins,
        z.basins$lmax_basins,
        z.basins$lmin_basins,
        if (is.null(vertex.mass)) NULL else as.numeric(vertex.mass),
        options,
        PACKAGE = "gflow"
    )

    ## Add metadata
    result$n_vertices <- length(y.hat)
    result$polarity_scale <- polarity.scale
    result$epsilon <- epsilon

    ## Ensure class is set (should already be set in C++)
    class(result) <- c("gfcor", "list")

    return(result)
}


#' Print Method for Gradient Flow Correlation
#'
#' @description
#' Displays a summary of gradient flow correlation results.
#'
#' @param x Object of class \code{"gfcor"} from \code{gfcor()}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.gfcor <- function(x, ...) {
    cat("Gradient Flow Correlation Analysis\n")
    cat("===================================\n\n")

    cat("Global Association Measures:\n")
    cat(sprintf("  Polarity concordance (A_pol):  %+.4f\n", x$global$A_pol))
    cat(sprintf("  Sign concordance (kappa_pol):  %+.4f\n", x$global$kappa_pol))
    cat("\n")

    cat("Vertex Counts:\n")
    cat(sprintf("  Positive association: %d\n", x$global$n_positive))
    cat(sprintf("  Negative association: %d\n", x$global$n_negative))
    cat(sprintf("  Zero association:     %d\n", x$global$n_zero))
    cat(sprintf("  Invalid (flat):       %d\n", x$global$n_invalid))
    cat("\n")

    cat("Basin Structure:\n")
    cat(sprintf("  Y maxima: %d, Y minima: %d\n",
                length(x$basin_character$chi_y_max),
                length(x$basin_character$chi_y_min)))
    cat(sprintf("  Z maxima: %d, Z minima: %d\n",
                length(x$basin_character$chi_z_max),
                length(x$basin_character$chi_z_min)))
    cat("\n")

    cat(sprintf("Settings: polarity.scale = '%s', epsilon = %.2e\n",
                x$polarity_scale, x$epsilon))

    invisible(x)
}


#' Summary Method for Gradient Flow Correlation
#'
#' @description
#' Generates a detailed summary of gradient flow correlation results including
#' basin-level statistics.
#'
#' @param object Object of class \code{"gfcor"} from \code{gfcor()}.
#' @param ... Additional arguments (unused).
#'
#' @return A list of class \code{"summary.gfcor"} containing summary tables.
#'
#' @export
summary.gfcor <- function(object, ...) {
    result <- list()

    ## Global summary
    result$global <- data.frame(
        measure = c("A_pol", "kappa_pol"),
        value = c(object$global$A_pol, object$global$kappa_pol),
        description = c("Polarity concordance", "Sign concordance")
    )

    ## Vertex distribution
    n.valid <- object$global$n_positive + object$global$n_negative + object$global$n_zero
    result$vertex_distribution <- data.frame(
        category = c("positive", "negative", "zero", "invalid"),
        count = c(object$global$n_positive, object$global$n_negative,
                  object$global$n_zero, object$global$n_invalid),
        proportion = c(object$global$n_positive, object$global$n_negative,
                       object$global$n_zero, object$global$n_invalid) / object$n_vertices
    )

    ## Basin character summary
    bc <- object$basin_character

    if (length(bc$chi_y_max) > 0) {
        result$y_max_basins <- data.frame(
            basin = seq_along(bc$chi_y_max),
            chi = bc$chi_y_max,
            mass = bc$mass_y_max
        )
    }

    if (length(bc$chi_y_min) > 0) {
        result$y_min_basins <- data.frame(
            basin = seq_along(bc$chi_y_min),
            chi = bc$chi_y_min,
            mass = bc$mass_y_min
        )
    }

    class(result) <- c("summary.gfcor", "list")
    return(result)
}


#' Print Method for Gradient Flow Correlation Summary
#'
#' @param x Object of class \code{"summary.gfcor"}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.gfcor <- function(x, ...) {
    cat("Summary: Gradient Flow Correlation Analysis\n")
    cat("============================================\n\n")

    cat("Global Measures:\n")
    print(x$global, row.names = FALSE)
    cat("\n")

    cat("Vertex Distribution:\n")
    print(x$vertex_distribution, row.names = FALSE)
    cat("\n")

    if (!is.null(x$y_max_basins)) {
        cat("Y Maximum Basins (association character chi):\n")
        print(x$y_max_basins, row.names = FALSE)
        cat("\n")
    }

    if (!is.null(x$y_min_basins)) {
        cat("Y Minimum Basins (association character chi):\n")
        print(x$y_min_basins, row.names = FALSE)
        cat("\n")
    }

    invisible(x)
}

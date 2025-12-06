## Local Slope (Asymmetric Association) Measures
##
## Functions for computing asymmetric local association measures between a
## directing function y and response function z defined on graph vertices.
## Unlike symmetric local correlation, these measures treat y as the directing
## variable and z as the response.

#' Gradient-Restricted Local Slope
#'
#' @description
#' Computes the local slope of z with respect to y along the gradient direction
#' of y at each vertex. This asymmetric measure captures "for each unit increase
#' in y along its steepest direction, how much does z change?"
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#'   Element i contains the neighbors of vertex i.
#' @param weight.list List of numeric vectors containing edge weights.
#'   Must have same structure as adj.list.
#' @param y Numeric vector of directing function values (length = number of vertices).
#' @param z Numeric vector of response function values (length = number of vertices).
#' @param type Character scalar specifying the type of measure:
#'   \itemize{
#'     \item "slope": Raw ratio Δz/Δy along gradient edge (unbounded)
#'     \item "normalized": Sigmoid-normalized ratio \eqn{\tanh(\alpha \Delta z/\Delta y)} (bounded to \eqn{[-1,1]})
#'     \item "sign": Sign of \eqn{\Delta z} along gradient direction \eqn{{-1, 0, +1}}
#'   }
#' @param y.diff.type Character scalar specifying edge difference type for y:
#'   \itemize{
#'     \item "difference": Standard differences (continuous data)
#'     \item "logratio": Log-ratios (compositional data)
#'   }
#' @param z.diff.type Character scalar specifying edge difference type for z.
#' @param epsilon Numeric scalar for pseudocount in log-ratios.
#'   If 0 (default), computed adaptively.
#' @param sigmoid.alpha Numeric scalar for sigmoid normalization scale.
#'   If 0 (default), auto-calibrated from median absolute slope.
#' @param ascending Logical. If TRUE (default), use ascending gradient
#'   (gradient toward higher y values). If FALSE, use descending gradient.
#' @param instrumented Logical. If TRUE, return full diagnostic output.
#'   If FALSE (default), return only coefficients.
#'
#' @return If instrumented=FALSE, returns a numeric vector of local slope
#'   coefficients. If instrumented=TRUE, returns a list with components:
#'   \describe{
#'     \item{vertex.coefficients}{Numeric vector of local slope at each vertex}
#'     \item{gradient.neighbors}{Integer vector of gradient edge neighbors (NA if extremum)}
#'     \item{gradient.delta.y}{Numeric vector of \eqn{\Delta y} along gradient edge}
#'     \item{gradient.delta.z}{Numeric vector of \eqn{\Delta z} along gradient edge}
#'     \item{is.local.extremum}{Logical vector indicating local extrema of y}
#'     \item{n.local.maxima}{Count of local maxima}
#'     \item{n.local.minima}{Count of local minima}
#'     \item{sigmoid.alpha}{Sigmoid scale used (for normalized type)}
#'     \item{mean.coefficient}{Mean of non-extremum coefficients}
#'     \item{median.coefficient}{Median of non-extremum coefficients}
#'   }
#'
#' @details
#' The gradient edge at vertex v is the edge from v to the neighbor u that
#' maximizes the difference \eqn{\Delta y = y(u) - y(v)} (for ascending=TRUE) or minimizes
#' it (for ascending=FALSE). At local extrema of y, no gradient edge exists
#' and the coefficient is set to 0.
#'
#' The three measure types capture different aspects:
#' \itemize{
#'   \item "slope": The raw regression coefficient \eqn{\Delta z/\Delta y}. Unbounded, can be
#'     sensitive near extrema where \Delta y is small.
#'   \item "normalized": Bounded to \eqn{[-1,1]} via sigmoid. More stable near extrema.
#'     The scale \eqn{\alpha} is auto-calibrated so median slopes map to ~0.5.
#'   \item "sign": Most robust but loses magnitude information. Answers only
#'     "does z increase or decrease along ∇y?"
#' }
#'
#' @section Relationship to Local Correlation:
#' Unlike symmetric local correlation lcor(y,z), the local slope is asymmetric:
#' lslope(z; y) ≠ lslope(y; z) in general. The slope measures the response of
#' z to changes in y, treating y as the independent variable.
#'
#' @section Application to Secondary Phylotypes:
#' For identifying secondary phylotypes associated with sPTB within Lactobacillus-
#' dominated samples, use the smoothed sPTB prevalence as y and phylotype
#' abundances as z. Positive slopes indicate phylotypes that increase as sPTB
#' risk increases along the gradient direction.
#'
#' @examples
#' \dontrun{
#' ## Simple example
#' adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
#' weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
#' y <- c(0.1, 0.5, 0.9)  # sPTB prevalence
#' z <- c(0.8, 0.4, 0.1)  # Secondary phylotype abundance
#'
#' ## Compute gradient-restricted slope
#' result <- lslope.gradient(adj.list, weight.list, y, z,
#'                           type = "normalized",
#'                           y.diff.type = "difference",
#'                           z.diff.type = "logratio",
#'                           instrumented = TRUE)
#'
#' ## Examine results
#' print(result$vertex.coefficients)
#' print(result$gradient.neighbors)
#' }
#'
#' @seealso \code{\link{lslope.neighborhood}} for neighborhood-based regression,
#'   \code{\link{lcor}} for symmetric local correlation
#'
#' @export
lslope.gradient <- function(adj.list,
                            weight.list,
                            y,
                            z,
                            type = c("normalized", "slope", "sign"),
                            y.diff.type = c("difference", "logratio"),
                            z.diff.type = c("difference", "logratio"),
                            epsilon = 0,
                            sigmoid.alpha = 0,
                            ascending = TRUE,
                            instrumented = FALSE) {

    ## Match arguments
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Validate inputs
    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }
    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }
    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }
    if (!is.numeric(z)) {
        stop("z must be a numeric vector")
    }
    if (!is.numeric(epsilon) || length(epsilon) != 1) {
        stop("epsilon must be a single numeric value")
    }
    if (!is.numeric(sigmoid.alpha) || length(sigmoid.alpha) != 1) {
        stop("sigmoid.alpha must be a single numeric value")
    }
    if (!is.logical(ascending) || length(ascending) != 1) {
        stop("ascending must be a single logical value")
    }

    n.vertices <- length(adj.list)

    if (length(y) != n.vertices) {
        stop("Length of y must equal number of vertices")
    }
    if (length(z) != n.vertices) {
        stop("Length of z must equal number of vertices")
    }
    if (length(weight.list) != n.vertices) {
        stop("Length of weight.list must equal number of vertices")
    }

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call C++ function
    if (instrumented) {
        result <- .Call(
            "S_lslope_gradient_instrumented",
            adj.list.0,
            weight.list,
            as.numeric(y),
            as.numeric(z),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(sigmoid.alpha),
            as.logical(ascending),
            PACKAGE = "gflow"
        )

        ## Add class for method dispatch
        class(result) <- c("lslope_gradient_result", "list")

        ## Add metadata
        attr(result, "type") <- type
        attr(result, "y.diff.type") <- y.diff.type
        attr(result, "z.diff.type") <- z.diff.type
        attr(result, "ascending") <- ascending

        return(result)
    } else {
        result <- .Call(
            "S_lslope_gradient",
            adj.list.0,
            weight.list,
            as.numeric(y),
            as.numeric(z),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(sigmoid.alpha),
            as.logical(ascending),
            PACKAGE = "gflow"
        )

        return(result)
    }
}


#' Neighborhood Local Regression Coefficient
#'
#' @description
#' Computes the local regression coefficient of z on y using all edges in each
#' vertex's neighborhood. This is the asymmetric analog of local correlation,
#' satisfying β_loc(z; y) = lcor(y, z) × sd_loc(z) / sd_loc(y).
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#' @param weight.list List of numeric vectors containing edge weights.
#' @param y Numeric vector of directing function values.
#' @param z Numeric vector of response function values.
#' @param weight.type Character scalar specifying weighting scheme:
#'   \itemize{
#'     \item "unit": Equal weights (w_e = 1)
#'     \item "derivative": Geometric weights (w_e = 1/length^2)
#'   }
#' @param y.diff.type Character scalar specifying edge difference type for y.
#' @param z.diff.type Character scalar specifying edge difference type for z.
#' @param epsilon Numeric scalar for pseudocount in log-ratios.
#'
#' @return A list with components:
#'   \describe{
#'     \item{vertex.coefficients}{Numeric vector of local regression coefficients}
#'     \item{sd.y}{Numeric vector of local standard deviation of y}
#'     \item{sd.z}{Numeric vector of local standard deviation of z}
#'     \item{lcor}{Numeric vector of local correlation (for reference)}
#'     \item{mean.coefficient}{Mean coefficient across all vertices}
#'     \item{median.coefficient}{Median coefficient across all vertices}
#'   }
#'
#' @details
#' The neighborhood local regression coefficient at vertex v is:
#' \deqn{\beta_{loc}(z; y, w)(v) = \frac{\sum w_e \Delta_e y \cdot \Delta_e z}{\sum w_e (\Delta_e y)^2}}
#'
#' Unlike the gradient-restricted slope which uses only the single gradient edge,
#' this measure uses all edges in the neighborhood, providing a more robust but
#' less direction-specific measure of local association.
#'
#' The relationship β = ρ × σ_z/σ_y holds locally: the regression coefficient
#' equals the correlation times the ratio of local standard deviations.
#'
#' @examples
#' \dontrun{
#' result <- lslope.neighborhood(adj.list, weight.list, y, z,
#'                               weight.type = "derivative",
#'                               y.diff.type = "difference",
#'                               z.diff.type = "logratio")
#'
#' ## Verify relationship with local correlation
#' beta.from.lcor <- result$lcor * result$sd.z / result$sd.y
#' all.equal(result$vertex.coefficients, beta.from.lcor)
#' }
#'
#' @seealso \code{\link{lslope.gradient}} for gradient-restricted slope,
#'   \code{\link{lcor}} for symmetric local correlation
#'
#' @export
lslope.neighborhood <- function(adj.list,
                                weight.list,
                                y,
                                z,
                                weight.type = c("derivative", "unit"),
                                y.diff.type = c("difference", "logratio"),
                                z.diff.type = c("difference", "logratio"),
                                epsilon = 0) {

    ## Match arguments
    weight.type <- match.arg(weight.type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Validate inputs
    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }
    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }
    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }
    if (!is.numeric(z)) {
        stop("z must be a numeric vector")
    }
    if (!is.numeric(epsilon) || length(epsilon) != 1) {
        stop("epsilon must be a single numeric value")
    }

    n.vertices <- length(adj.list)

    if (length(y) != n.vertices) {
        stop("Length of y must equal number of vertices")
    }
    if (length(z) != n.vertices) {
        stop("Length of z must equal number of vertices")
    }
    if (length(weight.list) != n.vertices) {
        stop("Length of weight.list must equal number of vertices")
    }

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call C++ function
    result <- .Call(
        "S_lslope_neighborhood",
        adj.list.0,
        weight.list,
        as.numeric(y),
        as.numeric(z),
        as.character(weight.type),
        as.character(y.diff.type),
        as.character(z.diff.type),
        as.numeric(epsilon),
        PACKAGE = "gflow"
    )

    ## Add class for method dispatch
    class(result) <- c("lslope_neighborhood_result", "list")

    ## Add metadata
    attr(result, "weight.type") <- weight.type
    attr(result, "y.diff.type") <- y.diff.type
    attr(result, "z.diff.type") <- z.diff.type
    attr(result, "n.vertices") <- n.vertices

    return(result)
}


#' Print Method for Gradient-Restricted Local Slope Results
#'
#' @param x An object of class "lslope_gradient_result"
#' @param ... Additional arguments (ignored)
#' @export
print.lslope_gradient_result <- function(x, ...) {
    cat("Gradient-Restricted Local Slope Result\n")
    cat("=======================================\n\n")

    cat("Parameters:\n")
    cat("  Type:", attr(x, "type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n")
    cat("  Ascending gradient:", attr(x, "ascending"), "\n\n")

    cat("Results:\n")
    cat("  Number of vertices:", x$n.vertices, "\n")
    cat("  Local maxima:", x$n.local.maxima, "\n")
    cat("  Local minima:", x$n.local.minima, "\n")
    cat("  Mean coefficient:", round(x$mean.coefficient, 4), "\n")
    cat("  Median coefficient:", round(x$median.coefficient, 4), "\n")

    if (attr(x, "type") == "normalized") {
        cat("  Sigmoid alpha:", round(x$sigmoid.alpha, 4), "\n")
    }

    cat("\nCoefficient distribution (non-extrema):\n")
    non.extrema <- x$vertex.coefficients[!x$is.local.extremum]
    if (length(non.extrema) > 0) {
        print(summary(non.extrema))
    }

    invisible(x)
}


#' Print Method for Neighborhood Local Slope Results
#'
#' @param x An object of class "lslope_neighborhood_result"
#' @param ... Additional arguments (ignored)
#' @export
print.lslope_neighborhood_result <- function(x, ...) {
    cat("Neighborhood Local Regression Coefficient Result\n")
    cat("=================================================\n\n")

    cat("Parameters:\n")
    cat("  Weight type:", attr(x, "weight.type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n\n")

    cat("Results:\n")
    cat("  Number of vertices:", attr(x, "n.vertices"), "\n")
    cat("  Mean coefficient:", round(x$mean.coefficient, 4), "\n")
    cat("  Median coefficient:", round(x$median.coefficient, 4), "\n")

    cat("\nCoefficient distribution:\n")
    print(summary(x$vertex.coefficients))

    cat("\nLocal correlation summary:\n")
    print(summary(x$lcor))

    invisible(x)
}

#' Local Correlation Coefficients with Flexible Edge Difference Types
#'
#' Compute vertex-level correlation coefficients measuring the alignment of
#' directional changes between two functions on a graph, with support for both
#' continuous (standard differences) and compositional (log-ratios) data.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#'   Element i contains the neighbors of vertex i.
#' @param weight.list List of numeric vectors containing edge weights.
#'   Must have same structure as adj.list.
#' @param y Numeric vector of response function values (length = number of vertices).
#' @param z Numeric vector of feature function values (length = number of vertices).
#' @param type Character scalar specifying weighting scheme:
#'   \itemize{
#'     \item "unit": Equal weights (w_e = 1)
#'     \item "derivative": Geometric weights (w_e = 1/length^2)
#'     \item "sign": Sign-based (typically not used with correlation)
#'   }
#' @param y.diff.type Character scalar specifying edge difference type for y:
#'   \itemize{
#'     \item "difference": Standard differences (continuous data)
#'     \item "logratio": Log-ratios (compositional data)
#'   }
#' @param z.diff.type Character scalar specifying edge difference type for z.
#'   Same options as y.diff.type.
#' @param epsilon Numeric scalar for pseudocount in log-ratios.
#'   If 0 (default), computed adaptively as 1e-6 * min(non-zero values).
#' @param winsorize.quantile Numeric scalar for winsorization.
#'   If 0 (default), no winsorization. If > 0 (e.g., 0.05), clips edge
#'   differences to [q, 1-q] percentiles for robustness.
#'
#' @return A list with components:
#'   \describe{
#'     \item{vertex.coefficients}{Numeric vector of correlation coefficients at each vertex}
#'     \item{vertex.delta.y}{List of numeric vectors, edge differences for y at each vertex}
#'     \item{vertex.delta.z}{List of numeric vectors, edge differences for z at each vertex}
#'     \item{vertex.weights}{List of numeric vectors, edge weights at each vertex}
#'     \item{all.delta.y}{Numeric vector of all y edge differences (pre-winsorization)}
#'     \item{all.delta.z}{Numeric vector of all z edge differences (pre-winsorization)}
#'     \item{y.lower}{Lower winsorization bound for y (-Inf if none)}
#'     \item{y.upper}{Upper winsorization bound for y (+Inf if none)}
#'     \item{z.lower}{Lower winsorization bound for z (-Inf if none)}
#'     \item{z.upper}{Upper winsorization bound for z (+Inf if none)}
#'   }
#'
#' @details
#' The local correlation coefficient at vertex v measures the alignment of
#' directional changes in y and z within v's neighborhood:
#'
#' \deqn{lcor(y,z)(v) = \frac{\sum w_e \Delta_e y \cdot \Delta_e z}
#'                           {\sqrt{\sum w_e (\Delta_e y)^2} \sqrt{\sum w_e (\Delta_e z)^2}}}
#'
#' The key innovation is flexible edge difference types:
#' \itemize{
#'   \item Standard differences: \eqn{\Delta_e f = f(u) - f(v)}
#'   \item Log-ratios: \eqn{\Delta_e f = \log((f(u) + \epsilon) / (f(v) + \epsilon))}
#' }
#'
#' This enables appropriate treatment of different data types:
#' \itemize{
#'   \item Continuous vs. continuous: Use "difference" for both
#'   \item Continuous vs. compositional: Use "difference" for continuous, "logratio" for compositional
#'   \item Compositional vs. compositional: Use "logratio" for both
#' }
#'
#' Log-ratios are appropriate for compositional data (relative abundances,
#' proportions) because they respect the simplex geometry via the Aitchison
#' distance. The correlation of log-ratio gradients measures directional
#' alignment in the natural geometry for compositional data.
#'
#' @examples
#' \dontrun{
#' # Continuous response vs. compositional feature
#' library(gflow)
#'
#' # Build graph (example: 3-vertex triangle)
#' adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
#' weight.list <- list(c(1, 1), c(1, 1), c(1, 1))
#'
#' # Data
#' severity <- c(0.1, 0.5, 0.9)      # Continuous (disease severity)
#' abundance <- c(0.01, 0.10, 0.50)  # Compositional (bacterial abundance)
#'
#' # Compute local correlation
#' result <- lcor(
#'   adj.list, weight.list,
#'   severity, abundance,
#'   type = "derivative",
#'   y.diff.type = "difference",   # Continuous
#'   z.diff.type = "logratio",     # Compositional
#'   epsilon = 0,                   # Adaptive
#'   winsorize.quantile = 0         # No winsorization
#' )
#'
#' # Examine results
#' print(result$vertex.coefficients)
#' hist(result$all.delta.z, main = "Distribution of log-ratios")
#'
#' # With winsorization for robustness
#' result.robust <- lcor(
#'   adj.list, weight.list,
#'   severity, abundance,
#'   type = "derivative",
#'   y.diff.type = "difference",
#'   z.diff.type = "logratio",
#'   epsilon = 1e-6,
#'   winsorize.quantile = 0.05     # Clip to 5th/95th percentiles
#' )
#' }
#'
#' @seealso \code{\link{comono}} for the standard difference-only version
#'
#' @export
lcor <- function(adj.list,
                 weight.list,
                 y,
                 z,
                 type = c("derivative", "unit", "sign"),
                 y.diff.type = c("difference", "logratio"),
                 z.diff.type = c("difference", "logratio"),
                 epsilon = 0,
                 winsorize.quantile = 0,
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
    if (!is.numeric(winsorize.quantile) || length(winsorize.quantile) != 1) {
        stop("winsorize.quantile must be a single numeric value")
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
    result <- NULL
    if (instrumented) {
        result <- .Call(
            "S_lcor_instrumented",
            adj.list.0,
            weight.list,
            as.numeric(y),
            as.numeric(z),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(winsorize.quantile),
            PACKAGE = "gflow"
        )
    } else {
        result <- .Call(
            "S_lcor",
            adj.list.0,
            weight.list,
            as.numeric(y),
            as.numeric(z),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(winsorize.quantile),
            PACKAGE = "gflow"
        )
    }

    ## Add class for method dispatch
    class(result) <- c("lcor_result", "list")

    ## Add metadata
    attr(result, "type") <- type
    attr(result, "y.diff.type") <- y.diff.type
    attr(result, "z.diff.type") <- z.diff.type
    attr(result, "epsilon") <- epsilon
    attr(result, "winsorize.quantile") <- winsorize.quantile
    attr(result, "n.vertices") <- n.vertices

    return(result)
}

#' Print method for lcor results
#'
#' @param x An object of class "lcor_result"
#' @param ... Additional arguments (ignored)
#' @export
print.lcor_result <- function(x, ...) {
  cat("Local Correlation Result\n")
  cat("========================\n\n")

  cat("Parameters:\n")
  cat("  Weighting type:", attr(x, "type"), "\n")
  cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
  cat("  Z difference type:", attr(x, "z.diff.type"), "\n")
  cat("  Epsilon:", attr(x, "epsilon"), "\n")
  cat("  Winsorize quantile:", attr(x, "winsorize.quantile"), "\n\n")

  cat("Results:\n")
  cat("  Number of vertices:", attr(x, "n.vertices"), "\n")
  cat("  Mean coefficient:", mean(x$vertex.coefficients), "\n")
  cat("  Median coefficient:", median(x$vertex.coefficients), "\n")
  cat("  Range:", paste(range(x$vertex.coefficients), collapse = " to "), "\n\n")

  # Show winsorization bounds if applicable
  if (is.finite(x$y.lower)) {
    cat("Winsorization bounds:\n")
    cat("  Y: [", x$y.lower, ", ", x$y.upper, "]\n", sep = "")
    cat("  Z: [", x$z.lower, ", ", x$z.upper, "]\n", sep = "")
  }

  invisible(x)
}

#' Summary method for lcor results
#'
#' @param object An object of class "lcor_result"
#' @param ... Additional arguments (ignored)
#' @export
summary.lcor_result <- function(object, ...) {
  cat("Local Correlation Summary\n")
  cat("=========================\n\n")

  cat("Vertex coefficients:\n")
  print(summary(object$vertex.coefficients))

  if (length(object$all.delta.y) > 0) {
    cat("\nY edge differences (n =", length(object$all.delta.y), "):\n")
    print(summary(object$all.delta.y))
  }

  if (length(object$all.delta.z) > 0) {
    cat("\nZ edge differences (n =", length(object$all.delta.z), "):\n")
    print(summary(object$all.delta.z))
  }

  invisible(object)
}

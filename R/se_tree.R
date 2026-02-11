#' Build Spurious Extrema Tree for a Single Spurious Extremum
#'
#' @description
#' Constructs a tree structure capturing the hierarchical relationships between
#' spurious extrema through their basins of attraction. The tree alternates
#' between minima and maxima, terminating at non-spurious extrema.
#'
#' @param adj.list Graph adjacency list (1-based indexing).
#' @param weight.list Edge weight list.
#' @param y Numeric vector of function values.
#' @param root.vertex Vertex index of the spurious extremum to use as root (1-based).
#' @param spurious.min Integer vector of spurious minimum vertices (1-based).
#' @param spurious.max Integer vector of spurious maximum vertices (1-based).
#' @param basins A gfc_basins object from compute.gfc.basins().
#' @param compute.hr.support Logical; if TRUE, compute the harmonic repair support region.
#' @param verbose Logical; print progress information.
#'
#' @return An object of class "se_tree" containing:
#'   \describe{
#'     \item{root.vertex}{Root vertex index (1-based)}
#'     \item{root.is.maximum}{Logical; TRUE if root is a maximum}
#'     \item{classification}{One of "both_types", "only_min", "only_max", "no_terminals"}
#'     \item{nodes}{Named list of node structures}
#'     \item{ns.min.terminals}{Non-spurious minimum terminal vertices}
#'     \item{ns.max.terminals}{Non-spurious maximum terminal vertices}
#'     \item{hr.support.vertices}{Vertices in HR support region (if computed)}
#'     \item{hr.support.boundary}{Boundary of HR support region (if computed)}
#'   }
#'
#' @details
#' The SE tree is built by breadth-first search starting from the root:
#' \enumerate{
#'   \item Find opposite-polarity extrema within the root's basin
#'   \item Add them as children, mark as discovered
#'   \item For each child, recursively find opposite-polarity extrema in its basin
#'   \item Terminate branches at non-spurious extrema
#' }
#'
#' The tree classification determines the HR support strategy:
#' \itemize{
#'   \item both_types: Use branches to both NSMIN and NSMAX terminals
#'   \item only_min/only_max: Use branch to nearest NS terminal
#'   \item no_terminals: Use all nodes (isolated SE cluster)
#' }
#'
#' @examples
#' \dontrun{
#' # First compute basins
#' basins <- compute.gfc.basins(adj.list, weight.list, fitted.values)
#'
#' # Identify spurious extrema (e.g., from hop index analysis)
#' spurious.min <- c(725, 809, 1806)
#' spurious.max <- c(718, 815)
#'
#' # Build SE tree from a spurious minimum
#' se.tree <- build.se.tree(
#'     adj.list, weight.list, fitted.values,
#'     root.vertex = 725,
#'     spurious.min = spurious.min,
#'     spurious.max = spurious.max,
#'     basins = basins,
#'     compute.hr.support = TRUE
#' )
#'
#' print(se.tree)
#' }
#'
#' @export
build.se.tree <- function(adj.list,
                           weight.list,
                           y,
                           root.vertex,
                           spurious.min,
                           spurious.max,
                           basins,
                           compute.hr.support = TRUE,
                           verbose = TRUE) {

    ## Input validation
    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }

    n.vertices <- length(adj.list)

    if (!is.numeric(root.vertex) || length(root.vertex) != 1) {
        stop("root.vertex must be a single integer")
    }
    root.vertex <- as.integer(root.vertex)
    if (root.vertex < 1 || root.vertex > n.vertices) {
        stop(sprintf("root.vertex must be between 1 and %d", n.vertices))
    }

    if (!inherits(basins, "gfc_basins")) {
        stop("basins must be a gfc_basins object from compute.gfc.basins()")
    }

    spurious.min <- as.integer(spurious.min)
    spurious.max <- as.integer(spurious.max)

    ## Convert adjacency list to 0-based
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) integer(0) else as.integer(x - 1)
    })

    ## Call C++
    result <- .Call(
        "S_build_se_tree",
        adj.list.0based,
        weight.list,
        as.numeric(y),
        root.vertex,
        spurious.min,
        spurious.max,
        basins,
        as.logical(compute.hr.support),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    return(result)
}


#' Print Method for se_tree Objects
#'
#' @param x An se_tree object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.se_tree <- function(x, ...) {
    cat("Spurious Extrema Tree\n")
    cat("=====================\n")
    cat(sprintf("Root: %s at vertex %d\n",
                ifelse(x$root.is.maximum, "maximum", "minimum"),
                x$root.vertex))
    cat(sprintf("Classification: %s\n", x$classification))
    cat(sprintf("Nodes: %d\n", length(x$nodes)))
    cat(sprintf("Non-spurious terminals: %d minima, %d maxima\n",
                length(x$ns.min.terminals),
                length(x$ns.max.terminals)))

    if (length(x$hr.support.vertices) > 0) {
        cat(sprintf("\nHR Support Region:\n"))
        cat(sprintf("  Interior vertices: %d\n", length(x$hr.support.vertices)))
        cat(sprintf("  Boundary vertices: %d\n", length(x$hr.support.boundary)))
    }

    ## Print tree structure
    cat("\nTree structure:\n")
    for (name in names(x$nodes)) {
        node <- x$nodes[[name]]
        indent <- ""
        if (!is.na(node$parent)) {
            indent <- "  "
        }
        status <- ifelse(node$is.spurious, "spurious", "NON-SPURIOUS")
        cat(sprintf("%s%s (vertex %d, %s, basin: %d vertices)\n",
                    indent, name, node$vertex, status,
                    length(node$basin.vertices)))
    }

    invisible(x)
}

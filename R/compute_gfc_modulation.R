#' Compute Gradient Flow Modulation Weights
#'
#' Computes modulation factors for gradient flow based on vertex density
#' and/or edge length distribution. Returns modulated edge weights along
#' with intermediate objects for inspection and debugging.
#'
#' @param adj.list Graph adjacency list (1-based indexing).
#' @param weight.list Edge weight list (edge lengths).
#' @param modulation Character string specifying modulation type.
#'   \code{"density"} for vertex density modulation (edges weighted by
#'   target vertex density), \code{"edgelen"} for edge length modulation
#'   (edges weighted by length distribution density), or
#'   \code{"density_edgelen"} for combined modulation (product of both).
#' @param density Optional numeric vector of pre-computed density values.
#'   If \code{NULL}, computed from nearest neighbor distances (raw d1 inverse).
#'   For smoothed density, pass result from \code{compute.smoothed.density()}.
#' @param edgelen.bandwidth Bandwidth for edge length KDE. If \code{"auto"}
#'   (default), uses Sheather-Jones method via \code{bw.SJ()}. Can also be
#'   a positive numeric value for manual bandwidth, or a character string
#'   for other bandwidth selectors: \code{"nrd0"}, \code{"nrd"}, \code{"ucv"},
#'   \code{"bcv"}, or \code{"SJ"}.
#' @param normalize.weights Logical; if \code{TRUE} (default), normalize
#'   modulated weights so that for each vertex the maximum outgoing weight
#'   equals the original maximum.
#' @param verbose Logical; print progress information.
#'
#' @return An object of class \code{"gfc_modulation"} containing:
#'   \item{modulated.weight.list}{List of modulated edge weights.}
#'   \item{modulation}{The modulation type used.}
#'   \item{density}{Vertex density values (\code{NULL} if not computed).}
#'   \item{density.raw}{Logical; \code{TRUE} if density was computed
#'     internally (raw d1 inverse).}
#'   \item{edgelen.weights}{Named list of edge length weights, where names
#'     are \code{"v-u"} for edge from v to u (\code{NULL} if not computed).}
#'   \item{edgelen.bandwidth}{Bandwidth used for edge length KDE.}
#'   \item{edgelen.kde}{The full density object from R's \code{density()}
#'     function, useful for plotting and inspection.}
#'   \item{edgelen.lengths}{Vector of all unique edge lengths.}
#'   \item{edgelen.densities}{KDE density at each edge length (before
#'     normalization).}
#'   \item{n.vertices}{Number of vertices.}
#'   \item{n.edges}{Number of unique edges.}
#'
#' @details
#' The modulation factors are applied to edge weights as follows.
#' For an edge from v to u with original weight w(v,u), the modulated
#' weight is: for \code{"density"}, w'(v,u) = w(v,u) * rho(u); for
#' \code{"edgelen"}, w'(v,u) = w(v,u) * dl(v,u); for \code{"density_edgelen"},
#' w'(v,u) = w(v,u) * rho(u) * dl(v,u). Here rho(u) is the density at
#' vertex u and dl(v,u) is the edge length density weight.
#'
#' @seealso \code{\link{compute.gfc.basins}}, \code{\link{compute.smoothed.density}}
#'
#' @export
compute.gfc.modulation <- function(adj.list,
                                   weight.list,
                                   modulation = c("density", "edgelen", "density_edgelen"),
                                   density = NULL,
                                   edgelen.bandwidth = "auto",
                                   normalize.weights = TRUE,
                                   verbose = TRUE) {

    ## ========================================================================
    ## Input validation
    ## ========================================================================

    modulation <- match.arg(modulation)

    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }

    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }

    n.vertices <- length(adj.list)

    if (length(weight.list) != n.vertices) {
        stop("adj.list and weight.list must have the same length")
    }

    ## ========================================================================
    ## Initialize result
    ## ========================================================================

    result <- list(
        modulated.weight.list = vector("list", n.vertices),
        modulation = modulation,
        density = NULL,
        density.raw = FALSE,
        edgelen.weights = NULL,
        edgelen.bandwidth = NA_real_,
        edgelen.kde = NULL,
        edgelen.lengths = NULL,
        edgelen.densities = NULL,
        n.vertices = n.vertices,
        n.edges = 0
    )

    ## ========================================================================
    ## Compute density if needed
    ## ========================================================================

    use.density <- modulation %in% c("density", "density_edgelen")

    if (use.density) {
        if (is.null(density)) {
            if (verbose) {
                cat("Computing raw vertex density (d1^-1)...\n")
            }
            density <- compute.vertex.density(adj.list, weight.list, normalize = TRUE)
            result$density.raw <- TRUE
        } else {
            if (length(density) != n.vertices) {
                stop("density must have length equal to number of vertices")
            }
            result$density.raw <- FALSE
        }
        result$density <- density

        if (verbose) {
            cat(sprintf("  Density range: [%.4f, %.4f], mean: %.4f\n",
                        min(density), max(density), mean(density)))
        }
    }

    ## ========================================================================
    ## Compute edge length weights if needed
    ## ========================================================================

    use.edgelen <- modulation %in% c("edgelen", "density_edgelen")

    if (use.edgelen) {
        if (verbose) {
            cat("Computing edge length distribution weights...\n")
        }

        ## Collect all unique edges and their lengths
        edge.list <- list()
        edge.lengths <- numeric(0)

        for (v in seq_len(n.vertices)) {
            nbrs <- adj.list[[v]]
            weights <- weight.list[[v]]

            for (j in seq_along(nbrs)) {
                u <- nbrs[j]
                if (u > v) {  ## Avoid double-counting
                    edge.key <- paste(v, u, sep = "-")
                    edge.list[[edge.key]] <- c(v, u)
                    edge.lengths <- c(edge.lengths, weights[j])
                }
            }
        }

        n.edges <- length(edge.lengths)
        result$n.edges <- n.edges
        result$edgelen.lengths <- edge.lengths

        if (n.edges == 0) {
            warning("No edges found in graph")
            result$edgelen.weights <- list()
        } else {
            ## Determine bandwidth
            if (is.character(edgelen.bandwidth)) {
                if (edgelen.bandwidth == "auto") {
                    bw <- "SJ"
                } else {
                    bw <- edgelen.bandwidth
                }
            } else if (is.numeric(edgelen.bandwidth) && edgelen.bandwidth > 0) {
                bw <- edgelen.bandwidth
            } else {
                stop("edgelen.bandwidth must be 'auto', a bandwidth method name, ",
                     "or a positive numeric value")
            }

            ## Compute KDE using R's density()
            kde <- tryCatch(
                density(edge.lengths, bw = bw, n = 512),
                error = function(e) {
                    if (verbose) {
                        cat(sprintf("  Warning: bw.%s failed, falling back to nrd0\n",
                                    as.character(bw)))
                    }
                    density(edge.lengths, bw = "nrd0", n = 512)
                }
            )

            result$edgelen.bandwidth <- kde$bw
            result$edgelen.kde <- kde  ## Store full density object for inspection

            if (verbose) {
                cat(sprintf("  Bandwidth method: %s, value: %.4f\n",
                            as.character(bw), kde$bw))
            }

            ## Interpolate KDE density at each observed edge length
            kde.at.edges <- approx(kde$x, kde$y, xout = edge.lengths, rule = 2)$y
            result$edgelen.densities <- kde.at.edges

            ## Normalize to [0, 1]
            max.density <- max(kde.at.edges)
            normalized.densities <- kde.at.edges / max.density

            if (verbose) {
                cat(sprintf("  Edge lengths: [%.4f, %.4f], max density: %.4f\n",
                            min(edge.lengths), max(edge.lengths), max.density))
            }

            ## Build edge weight lookup (symmetric)
            edgelen.weights <- list()
            edge.keys <- names(edge.list)
            for (i in seq_along(edge.keys)) {
                key <- edge.keys[i]
                v <- edge.list[[key]][1]
                u <- edge.list[[key]][2]
                w <- normalized.densities[i]

                ## Store both directions
                edgelen.weights[[paste(v, u, sep = "-")]] <- w
                edgelen.weights[[paste(u, v, sep = "-")]] <- w
            }

            result$edgelen.weights <- edgelen.weights
        }
    }

    ## ========================================================================
    ## Compute modulated weights
    ## ========================================================================

    if (verbose) {
        cat("Computing modulated edge weights...\n")
    }

    for (v in seq_len(n.vertices)) {
        nbrs <- adj.list[[v]]
        orig.weights <- weight.list[[v]]
        n.nbrs <- length(nbrs)

        if (n.nbrs == 0) {
            result$modulated.weight.list[[v]] <- numeric(0)
            next
        }

        mod.weights <- numeric(n.nbrs)

        for (j in seq_len(n.nbrs)) {
            u <- nbrs[j]
            w <- orig.weights[j]

            ## Compute modulation factor
            mod.factor <- 1.0

            if (use.density) {
                mod.factor <- mod.factor * density[u]
            }

            if (use.edgelen) {
                edge.key <- paste(v, u, sep = "-")
                if (edge.key %in% names(result$edgelen.weights)) {
                    mod.factor <- mod.factor * result$edgelen.weights[[edge.key]]
                }
            }

            mod.weights[j] <- mod.factor
        }

        ## Normalize if requested (preserve max weight per vertex)
        if (normalize.weights && max(mod.weights) > 0) {
            ## Scale so max modulated weight equals max original weight
            ## This preserves the relative importance while applying modulation
            scale.factor <- max(orig.weights) / max(mod.weights)
            mod.weights <- mod.weights * scale.factor
        }

        result$modulated.weight.list[[v]] <- mod.weights
    }

    if (verbose) {
        cat("Modulation complete.\n")
    }

    class(result) <- "gfc_modulation"

    return(result)
}


#' Print Method for gfc_modulation Objects
#'
#' @param x A gfc_modulation object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.gfc_modulation <- function(x, ...) {

    cat("Gradient Flow Modulation\n")
    cat("========================\n")
    cat(sprintf("Type: %s\n", x$modulation))
    cat(sprintf("Vertices: %d, Edges: %d\n", x$n.vertices, x$n.edges))

    if (!is.null(x$density)) {
        cat(sprintf("\nDensity:\n"))
        cat(sprintf("  Source: %s\n",
                    ifelse(x$density.raw, "raw (d1^-1)", "user-provided")))
        cat(sprintf("  Range: [%.4f, %.4f]\n", min(x$density), max(x$density)))
        cat(sprintf("  Mean: %.4f, SD: %.4f\n", mean(x$density), sd(x$density)))
    }

    if (!is.null(x$edgelen.weights)) {
        cat(sprintf("\nEdge Length Distribution:\n"))
        cat(sprintf("  Bandwidth: %.4f\n", x$edgelen.bandwidth))
        cat(sprintf("  Length range: [%.4f, %.4f]\n",
                    min(x$edgelen.lengths), max(x$edgelen.lengths)))

        ## Find mode
        mode.idx <- which.max(x$edgelen.densities)
        cat(sprintf("  Mode: %.4f (length), %.4f (density)\n",
                    x$edgelen.lengths[mode.idx],
                    x$edgelen.densities[mode.idx]))
    }

    invisible(x)
}


#' Plot Diagnostics for Gradient Flow Modulation
#'
#' @param x A gfc_modulation object
#' @param type Type of plot: "density", "edgelen", or "both"
#' @param ... Additional arguments passed to plot functions
#'
#' @export
plot.gfc_modulation <- function(x,
                                   type = c("both", "density", "edgelen"),
                                   ...) {

    type <- match.arg(type)

    has.density <- !is.null(x$density)
    has.edgelen <- !is.null(x$edgelen.weights)

    if (type == "both" && has.density && has.edgelen) {
        par(mfrow = c(2, 2))
    } else if (type == "both" && (has.density || has.edgelen)) {
        par(mfrow = c(1, 2))
    }

    ## Density plots
    if ((type %in% c("both", "density")) && has.density) {
        ## Histogram
        hist(x$density, breaks = 30, col = "lightblue",
             main = "Vertex Density Distribution",
             xlab = expression(rho(v)),
             ...)

        ## Sorted density
        plot(sort(x$density), type = "l", col = "blue", lwd = 2,
             main = "Sorted Vertex Density",
             xlab = "Vertex rank", ylab = expression(rho(v)),
             ...)
        grid()
    }

    ## Edge length plots
    if ((type %in% c("both", "edgelen")) && has.edgelen) {
        ## Histogram of edge lengths with KDE overlay
        hist(x$edgelen.lengths, breaks = 30, col = "lightgreen",
             freq = FALSE,
             main = "Edge Length Distribution",
             xlab = "Edge length",
             ...)

        ## Overlay KDE curve
        if (!is.null(x$edgelen.kde)) {
            lines(x$edgelen.kde, col = "darkgreen", lwd = 2)
        }

        ## Normalized weights vs edge length
        ord <- order(x$edgelen.lengths)
        plot(x$edgelen.lengths[ord],
             x$edgelen.densities[ord] / max(x$edgelen.densities),
             type = "p", pch = 16, cex = 0.5, col = "darkgreen",
             main = "Edge Length Modulation Weights",
             xlab = "Edge length", ylab = "Normalized weight (dl)",
             ...)
        grid()

        ## Add KDE-based curve
        if (!is.null(x$edgelen.kde)) {
            kde.norm <- x$edgelen.kde$y / max(x$edgelen.kde$y)
            lines(x$edgelen.kde$x, kde.norm, col = "red", lwd = 2)
        }
    }

    if (type == "both") {
        par(mfrow = c(1, 1))
    }

    invisible(x)
}


#' Extract Modulation Factor for a Specific Edge
#'
#' @param mod A gfc_modulation object
#' @param v Source vertex (1-based)
#' @param u Target vertex (1-based)
#'
#' @return Named list with density, edgelen, and combined modulation factors
#'
#' @export
edge.modulation <- function(mod, v, u) {

    if (!inherits(mod, "gfc_modulation")) {
        stop("mod must be a gfc_modulation object")
    }

    result <- list(
        v = v,
        u = u,
        density.factor = NA_real_,
        edgelen.factor = NA_real_,
        combined.factor = NA_real_
    )

    ## Density factor (at target vertex)
    if (!is.null(mod$density)) {
        result$density.factor <- mod$density[u]
    }

    ## Edge length factor
    if (!is.null(mod$edgelen.weights)) {
        edge.key <- paste(v, u, sep = "-")
        if (edge.key %in% names(mod$edgelen.weights)) {
            result$edgelen.factor <- mod$edgelen.weights[[edge.key]]
        }
    }

    ## Combined factor
    combined <- 1.0
    if (!is.na(result$density.factor)) {
        combined <- combined * result$density.factor
    }
    if (!is.na(result$edgelen.factor)) {
        combined <- combined * result$edgelen.factor
    }
    result$combined.factor <- combined

    return(result)
}

#' Compute and Refine Basins of Attraction
#'
#' This function computes the initial basins of attraction for a gradient flow
#' and applies a sequence of filtering and merging operations to produce a
#' refined basin structure. The refinement process includes filtering by
#' relative values, clustering and merging similar extrema, and removing
#' extrema based on geometric characteristics.
#'
#' @param adj.list Adjacency list representation of the graph structure
#' @param edge.length.list List of edge lengths corresponding to the adjacency list
#' @param fitted.values Numeric vector of fitted values at each vertex
#' @param min.rel.value.max Minimum relative value threshold for maxima (default: 1.1)
#' @param max.rel.value.min Maximum relative value threshold for minima (default: 0.9)
#' @param max.overlap.threshold Overlap threshold for clustering maxima (default: 0.15)
#' @param min.overlap.threshold Overlap threshold for clustering minima (default: 0.15)
#' @param p.mean.hopk.dist.threshold Threshold for mean hop-k distance percentile (default: 0.9)
#' @param p.deg.threshold Threshold for degree percentile (default: 0.9)
#' @param hop.k Parameter for hop-k distance calculation in summary (default: 2)
#' @param apply.relvalue.filter Logical indicating whether to apply relative value filtering (default: TRUE)
#' @param apply.maxima.clustering Logical indicating whether to cluster and merge maxima (default: TRUE)
#' @param apply.minima.clustering Logical indicating whether to cluster and merge minima (default: TRUE)
#' @param apply.geometric.filter Logical indicating whether to apply geometric filtering (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: FALSE)
#'
#' @return A list with two components:
#'   \item{basins}{The refined basins object containing:
#'     \describe{
#'       \item{lmax_basins}{Named list of maximum basins}
#'       \item{lmin_basins}{Named list of minimum basins}
#'       \item{y}{Vector of function values}
#'       \item{adj.list}{Adjacency list (graph structure)}
#'       \item{edge.length.list}{Edge lengths (if provided)}
#'       \item{n.vertices}{Number of vertices}
#'     }
#'     The basin lists are named with extremum labels from the final summary,
#'     enabling direct access and simplifying merge operations.
#'   }
#'   \item{summary}{A data frame summarizing characteristics of refined basins}
#'
#' @details
#' The refinement process consists of several stages. First, the function computes
#' the initial basins of attraction using the gradient flow structure. These basins
#' capture the regions of influence for each local extremum in the fitted surface.
#'
#' The relative value filtering stage removes extrema whose values are too close
#' to the median of the fitted values. Maxima with relative values below
#' \code{min.rel.value.max} and minima with relative values above
#' \code{max.rel.value.min} are eliminated. This step focuses subsequent analysis
#' on prominent features of the fitted surface.
#'
#' The clustering stage identifies groups of similar extrema based on basin
#' overlap. When basins share a substantial portion of their vertices, their
#' corresponding extrema likely represent the same underlying feature. The
#' \code{max.overlap.threshold} and \code{min.overlap.threshold} parameters
#' control the sensitivity of this clustering for maxima and minima respectively.
#'
#' After clustering, the merging stage combines basins within each cluster into
#' a single basin. The merged basin inherits the extremum with the most extreme
#' value within the cluster. This reduces redundancy in the basin structure while
#' preserving the essential geometric features.
#'
#' The geometric filtering stage removes extrema whose basins exhibit unusual
#' structural characteristics. Basins with high values of mean hop-k distance or
#' degree percentiles may represent spurious features or boundary artifacts. The
#' thresholds \code{p.mean.hopk.dist.threshold} and \code{p.deg.threshold}
#' determine which basins to retain based on these geometric measures.
#'
#' Each filtering stage can be disabled by setting the corresponding logical
#' parameter to FALSE, allowing users to customize the refinement pipeline
#' for specific applications.
#'
#' @examples
#' \dontrun{
#' # Compute refined basins with default parameters
#' result <- compute.refined.basins(adj.list, edge.length.list, 
#'                                  fitted.values, verbose = TRUE)
#'
#' # Access the refined basins and summary
#' refined.basins <- result$basins
#' basin.summary <- result$summary
#'
#' # Use custom thresholds for more aggressive filtering
#' result <- compute.refined.basins(adj.list, edge.length.list, fitted.values,
#'                                  min.rel.value.max = 1.2,
#'                                  max.rel.value.min = 0.8,
#'                                  p.mean.hopk.dist.threshold = 0.85,
#'                                  verbose = TRUE)
#'
#' # Skip clustering steps
#' result <- compute.refined.basins(adj.list, edge.length.list, fitted.values,
#'                                  apply.maxima.clustering = FALSE,
#'                                  apply.minima.clustering = FALSE,
#'                                  verbose = TRUE)
#' }
#'
#' @export
compute.refined.basins <- function(adj.list,
                                   edge.length.list,
                                   fitted.values,
                                   min.rel.value.max = 1.1,
                                   max.rel.value.min = 0.9,
                                   max.overlap.threshold = 0.15,
                                   min.overlap.threshold = 0.15,
                                   p.mean.hopk.dist.threshold = 0.9,
                                   p.deg.threshold = 0.9,
                                   hop.k = 2,
                                   apply.relvalue.filter = TRUE,
                                   apply.maxima.clustering = TRUE,
                                   apply.minima.clustering = TRUE,
                                   apply.geometric.filter = TRUE,
                                   verbose = FALSE) {
  
    ## Validate inputs
    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }
    if (!is.list(edge.length.list)) {
        stop("edge.length.list must be a list")
    }
    if (!is.numeric(fitted.values)) {
        stop("fitted.values must be numeric")
    }
    if (length(adj.list) != length(edge.length.list)) {
        stop("adj.list and edge.length.list must have the same length")
    }
    if (length(fitted.values) != length(adj.list)) {
        stop("fitted.values must have the same length as adj.list")
    }

    ## Step 1: Compute initial basins of attraction
    if (verbose) {
        cat("Step 1: Computing initial basins of attraction...\n")
    }

    current.basins <- compute.basins.of.attraction(adj.list,
                                                   edge.length.list,
                                                   fitted.values)

    if (verbose) {
        initial.summary <- summary(current.basins, adj.list, edge.length.list)
        n.max <- sum(initial.summary$type == "max")
        n.min <- sum(initial.summary$type == "min")
        cat(sprintf("  Found %d maxima and %d minima\n", n.max, n.min))
    }

    ## Step 2: Filter by relative values
    if (apply.relvalue.filter) {
        if (verbose) {
            cat(sprintf("\nStep 2: Filtering by relative values (max >= %.2f, min <= %.2f)...\n",
                        min.rel.value.max, max.rel.value.min))
        }

        current.basins <- filter.basins.by.relvalue(current.basins,
                                                    min.rel.value.max = min.rel.value.max,
                                                    max.rel.value.min = max.rel.value.min)

        if (verbose) {
            relvalue.summary <- summary(current.basins, adj.list, edge.length.list)
            n.max <- sum(relvalue.summary$type == "max")
            n.min <- sum(relvalue.summary$type == "min")
            cat(sprintf("  Retained %d maxima and %d minima\n", n.max, n.min))
        }
    } else {
        if (verbose) {
            cat("\nStep 2: Skipping relative value filtering\n")
        }
    }

    ## Step 3: Cluster and merge maxima
    if (apply.maxima.clustering) {
        if (verbose) {
            cat(sprintf("\nStep 3: Clustering maxima (overlap threshold = %.2f)...\n",
                        max.overlap.threshold))
        }

        max.clusters <- cluster.local.extrema(current.basins,
                                              adj.list,
                                              edge.length.list,
                                              extrema.type = "max",
                                              overlap.threshold = max.overlap.threshold)

        if (verbose) {
            n.clusters <- length(unique(max.clusters$cluster))
            n.singletons <- sum(table(max.clusters$cluster) == 1)
            n.merged <- n.clusters - n.singletons
            cat(sprintf("  Found %d clusters (%d singletons, %d to be merged)\n",
                        n.clusters, n.singletons, n.merged))
            cat("  Merging clustered maxima...\n")
        }

        current.basins <- merge.clustered.basins(current.basins,
                                                 max.clusters,
                                                 extrema.type = "max")

        if (verbose) {
            merged.max.summary <- summary(current.basins, adj.list, edge.length.list)
            n.max <- sum(merged.max.summary$type == "max")
            cat(sprintf("  Result: %d maxima after merging\n", n.max))
        }
    } else {
        if (verbose) {
            cat("\nStep 3: Skipping maxima clustering\n")
        }
    }

    ## Step 4: Cluster and merge minima
    if (apply.minima.clustering) {
        if (verbose) {
            cat(sprintf("\nStep 4: Clustering minima (overlap threshold = %.2f)...\n",
                        min.overlap.threshold))
        }

        min.clusters <- cluster.local.extrema(current.basins,
                                              adj.list,
                                              edge.length.list,
                                              extrema.type = "min",
                                              overlap.threshold = min.overlap.threshold)

        if (verbose) {
            n.clusters <- length(unique(min.clusters$cluster))
            n.singletons <- sum(table(min.clusters$cluster) == 1)
            n.merged <- n.clusters - n.singletons
            cat(sprintf("  Found %d clusters (%d singletons, %d to be merged)\n",
                        n.clusters, n.singletons, n.merged))
            cat("  Merging clustered minima...\n")
        }

        current.basins <- merge.clustered.basins(current.basins,
                                                 min.clusters,
                                                 extrema.type = "min")

        if (verbose) {
            merged.min.summary <- summary(current.basins, adj.list, edge.length.list)
            n.min <- sum(merged.min.summary$type == "min")
            cat(sprintf("  Result: %d minima after merging\n", n.min))
        }
    } else {
        if (verbose) {
            cat("\nStep 4: Skipping minima clustering\n")
        }
    }

    ## Step 5: Filter by geometric characteristics
    if (apply.geometric.filter) {
        if (verbose) {
            cat(sprintf("\nStep 5: Filtering by geometric characteristics (p.mean.hopk.dist < %.2f, p.deg < %.2f)...\n",
                        p.mean.hopk.dist.threshold, p.deg.threshold))
        }

        ## Get summary with geometric characteristics
        geom.summary <- summary(current.basins, adj.list, edge.length.list, hop.k = hop.k)

        ## Filter maxima
        max.basin.df <- geom.summary[geom.summary$type == "max", ]
        good.max.vertices <- max.basin.df$vertex[
                                              max.basin.df$p.mean.hopk.dist < p.mean.hopk.dist.threshold &
                                              max.basin.df$p.deg < p.deg.threshold
                                          ]

        if (verbose) {
            n.max.before <- nrow(max.basin.df)
            n.max.after <- length(good.max.vertices)
            cat(sprintf("  Maxima: %d -> %d (removed %d)\n",
                        n.max.before, n.max.after, n.max.before - n.max.after))
        }

        ## Filter minima
        min.basin.df <- geom.summary[geom.summary$type == "min", ]
        good.min.vertices <- min.basin.df$vertex[
                                              min.basin.df$p.mean.hopk.dist < p.mean.hopk.dist.threshold &
                                              min.basin.df$p.deg < p.deg.threshold
                                          ]

        if (verbose) {
            n.min.before <- nrow(min.basin.df)
            n.min.after <- length(good.min.vertices)
            cat(sprintf("  Minima: %d -> %d (removed %d)\n",
                        n.min.before, n.min.after, n.min.before - n.min.after))
        }

        ## Combine good extrema and filter basins
        good.extrema <- c(good.min.vertices, good.max.vertices)
        current.basins <- filter.basins(current.basins, good.extrema)

    } else {
        if (verbose) {
            cat("\nStep 5: Skipping geometric filtering\n")
        }
    }

    ## Generate final summary
    if (verbose) {
        cat(sprintf("\nGenerating final summary (hop.k = %d)...\n", hop.k))
    }

    final.summary <- summary(current.basins, adj.list, edge.length.list, hop.k = hop.k)

    ## Populate names in lmin_basins and lmax_basins using extrema labels
    if (verbose) {
        cat("Populating extrema labels in basins object...\n")
    }

    ## Extract vertex IDs for maxima and minima
    max.vertices <- final.summary$label[final.summary$type == "max"]
    min.vertices <- final.summary$label[final.summary$type == "min"]

    ## Name the basin list elements using vertex IDs
    if (length(max.vertices) > 0) {
        names(current.basins$lmax_basins) <- as.character(max.vertices)
    }

    if (length(min.vertices) > 0) {
        names(current.basins$lmin_basins) <- as.character(min.vertices)
    }

    ## Store graph structure in basins object for downstream operations
    ## This makes the basins object self-contained and eliminates the need
    ## to pass adj.list and edge.length.list to merge.two.extrema()
    if (verbose) {
        cat("Storing graph structure in basins object...\n")
    }

    ## Always store adj.list (required for merge operations)
    if (is.null(current.basins$adj.list)) {
        current.basins$adj.list <- adj.list
    }

    if (is.null(current.basins$edge.length.list)) {
        current.basins$edge.length.list <- edge.length.list
    }

    current.basins$hop.k <- hop.k
    current.basins$n.vertices <- current.basins$n_vertices
    current.basins$n_vertices <- NULL

    if (verbose) {
        n.max.final <- sum(final.summary$type == "max")
        n.min.final <- sum(final.summary$type == "min")
        cat("\nRefinement complete!\n")
        cat(sprintf("Final structure: %d maxima and %d minima\n", n.max.final, n.min.final))
    }

    ## Return both basins and summary
    result <- list(
        basins = current.basins,
        summary = final.summary
    )

    class(result) <- "basins_of_attraction"

    return(result)
}


#' Print Method for Refined Basins
#'
#' @param x A refined.basins object
#' @param ... Additional arguments passed to print methods
#'
#' @export
print.refined.basins <- function(x, ...) {
  cat("Refined Basins of Attraction\n")
  cat("============================\n\n")
  
  n.max <- sum(x$summary$type == "max")
  n.min <- sum(x$summary$type == "min")
  
  cat(sprintf("Number of maxima: %d\n", n.max))
  cat(sprintf("Number of minima: %d\n", n.min))
  cat(sprintf("Total extrema: %d\n\n", n.max + n.min))
  
  cat("Basin summary (first 10 rows):\n")
  print(head(x$summary, 10))
  
  if (nrow(x$summary) > 10) {
    cat(sprintf("\n... and %d more rows\n", nrow(x$summary) - 10))
  }
  
  invisible(x)
}

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
#' @param edge.length.quantile.thld Numeric value in (0, 1] specifying the quantile
#'   of the graph's edge length distribution to use as a threshold for basin
#'   construction. This parameter prevents "basin jumping" pathology by excluding
#'   long edges that may cause trajectories to skip over intermediate critical points.
#' @param min.rel.value.max Minimum relative value threshold for maxima (default: 1.1)
#' @param max.rel.value.min Maximum relative value threshold for minima (default: 0.9)
#' @param max.overlap.threshold Overlap threshold for clustering maxima (default: 0.15)
#' @param min.overlap.threshold Overlap threshold for clustering minima (default: 0.15)
#' @param p.mean.nbrs.dist.threshold Threshold for mean nearest neighbors distance percentile (default: 0.9)
#' @param p.mean.hopk.dist.threshold Threshold for mean hop-k distance percentile (default: 0.9)
#' @param p.deg.threshold Threshold for degree percentile (default: 0.9)
#' @param min.basin.size Minimum basin size threshold. Extrema with basins
#'   containing fewer than this many vertices are removed. Default is 1 (no filtering).
#' @param hop.k Parameter for hop-k distance calculation in summary (default: 2)
#' @param apply.relvalue.filter Logical indicating whether to apply relative value filtering (default: TRUE)
#' @param apply.maxima.clustering Logical indicating whether to cluster and merge maxima (default: TRUE)
#' @param apply.minima.clustering Logical indicating whether to cluster and merge minima (default: TRUE)
#' @param apply.geometric.filter Logical indicating whether to apply geometric filtering (default: TRUE)
#' @param with.trajectories Set to TRUE for the function to return gradient trajectories.
#' @param verbose Logical indicating whether to print progress messages (default: FALSE)
#'
#' @return A list of class \code{"basins_of_attraction"} with the following components:
#'   \describe{
#'     \item{basins}{The refined basins object containing:
#'       \describe{
#'         \item{lmax_basins}{Named list of maximum basins, where names correspond
#'           to labels (M1, M2, ...) from the summary. Each basin contains a
#'           \code{vertex} (the extremum location) and \code{basin_df} (data frame
#'           of basin vertices with gradient flow information).}
#'         \item{lmin_basins}{Named list of minimum basins, with analogous structure
#'           and labels (m1, m2, ...).}
#'         \item{y}{Numeric vector of function values at each vertex.}
#'         \item{adj.list}{Adjacency list representation of the graph structure.}
#'         \item{edge.length.list}{List of edge lengths corresponding to the
#'           adjacency list.}
#'         \item{n.vertices}{Integer giving the number of vertices in the graph.}
#'         \item{hop.k}{The hop distance parameter used for summary statistics.}
#'       }
#'     }
#'     \item{summary}{A data frame with one row per retained extremum, containing:
#'       label, vertex index, function value, relative value, extremum type,
#'       hop index, basin size, distance percentiles, degree, and degree percentile.
#'       Minima are labeled m1, m2, ... in order of increasing value; maxima are
#'       labeled M1, M2, ... in order of decreasing value.}
#'     \item{max.overlap.dist}{Symmetric matrix of pairwise overlap distances
#'       between maximum basins, with row and column names matching the summary
#'       labels. Returns \code{NULL} if fewer than two maxima remain after
#'       refinement.}
#'     \item{min.overlap.dist}{Symmetric matrix of pairwise overlap distances
#'       between minimum basins, with row and column names matching the summary
#'       labels. Returns \code{NULL} if fewer than two minima remain after
#'       refinement.}
#'   }
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
#' result <- compute.refined.basins(adj.list,
#'                                  edge.length.list,
#'                                  fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  verbose = TRUE)
#'
#' # Access the refined basins and summary
#' refined.basins <- result$basins
#' basin.summary <- result$summary
#'
#' # Use custom thresholds for more aggressive filtering
#' result <- compute.refined.basins(adj.list, edge.length.list, fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  min.rel.value.max = 1.2,
#'                                  max.rel.value.min = 0.8,
#'                                  p.mean.hopk.dist.threshold = 0.85,
#'                                  verbose = TRUE)
#'
#' # Skip clustering steps
#' result <- compute.refined.basins(adj.list, edge.length.list, fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  apply.maxima.clustering = FALSE,
#'                                  apply.minima.clustering = FALSE,
#'                                  verbose = TRUE)
#' }
#'
#' @export
compute.refined.basins <- function(adj.list,
                                   edge.length.list,
                                   fitted.values,
                                   edge.length.quantile.thld = 0.9,
                                   min.rel.value.max = 1.1,
                                   max.rel.value.min = 0.9,
                                   max.overlap.threshold = 0.15,
                                   min.overlap.threshold = 0.15,
                                   p.mean.nbrs.dist.threshold = 0.9,
                                   p.mean.hopk.dist.threshold = 0.9,
                                   p.deg.threshold = 0.9,
                                   min.basin.size = 1,
                                   hop.k = 2,
                                   apply.relvalue.filter = TRUE,
                                   apply.maxima.clustering = TRUE,
                                   apply.minima.clustering = TRUE,
                                   apply.geometric.filter = TRUE,
                                   with.trajectories = FALSE,
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
                                                   fitted.values,
                                                   edge.length.quantile.thld,
                                                   with.trajectories)

    if (verbose) {
        initial.summary <- summary(current.basins, adj.list, edge.length.list, hop.k)
        n.max <- sum(initial.summary$type == "max")
        n.min <- sum(initial.summary$type == "min")
        cat(sprintf("  Found %d maxima and %d minima\n", n.max, n.min))
        cat("initial.summary:\n")
        print(initial.summary)
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
            relvalue.summary <- summary(current.basins, adj.list, edge.length.list, hop.k)
            n.max <- sum(relvalue.summary$type == "max")
            n.min <- sum(relvalue.summary$type == "min")
            cat(sprintf("  Retained %d maxima and %d minima\n", n.max, n.min))
            cat("relvalue.summary:\n")
            print(relvalue.summary)
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
            merged.max.summary <- summary(current.basins, adj.list, edge.length.list, hop.k)
            n.max <- sum(merged.max.summary$type == "max")
            cat(sprintf("  Result: %d maxima after merging\n", n.max))
            cat("merged.max.summary:\n")
            print(merged.max.summary)
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
            merged.min.summary <- summary(current.basins, adj.list, edge.length.list, hop.k)
            n.min <- sum(merged.min.summary$type == "min")
            cat(sprintf("  Result: %d minima after merging\n", n.min))
            cat("merged.min.summary:\n")
            print(merged.min.summary)
        }
    } else {
        if (verbose) {
            cat("\nStep 4: Skipping minima clustering\n")
        }
    }

    ## Step 5: Filter by geometric characteristics and basin size
    if (apply.geometric.filter) {
        if (verbose) {
            cat(sprintf("\nStep 5: Filtering by geometric characteristics and basin size:\n"))
            cat(sprintf("  p.mean.nbrs.dist < %.2f, p.mean.hopk.dist < %.2f, p.deg < %.2f, basin.size >= %d ...\n",
                        p.mean.nbrs.dist.threshold, p.mean.hopk.dist.threshold,
                        p.deg.threshold, min.basin.size))
        }

        ## Get summary with geometric characteristics
        geom.summary <- summary(current.basins, adj.list, edge.length.list, hop.k = hop.k)

        ## Filter maxima
        max.basin.df <- geom.summary[geom.summary$type == "max", ]
        good.max.vertices <- max.basin.df$vertex[
                                              max.basin.df$p.mean.nbrs.dist < p.mean.nbrs.dist.threshold &
                                              max.basin.df$p.mean.hopk.dist < p.mean.hopk.dist.threshold &
                                              max.basin.df$p.deg < p.deg.threshold &
                                              max.basin.df$basin.size >= min.basin.size
                                          ]

        if (verbose) {
            n.max.before <- nrow(max.basin.df)
            n.max.after <- length(good.max.vertices)
            cat(sprintf("  Maxima: %d -> %d (removed %d)\n",
                        n.max.before, n.max.after, n.max.before - n.max.after))
            cat("geom.summary:\n")
            print(geom.summary)
        }

        ## Filter minima
        min.basin.df <- geom.summary[geom.summary$type == "min", ]
        good.min.vertices <- min.basin.df$vertex[
                                              min.basin.df$p.mean.hopk.dist < p.mean.hopk.dist.threshold &
                                              min.basin.df$p.deg < p.deg.threshold &
                                              min.basin.df$basin.size >= min.basin.size
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
        cat("final.summary:\n")
        print(final.summary)
    }

    ## Create lookup tables: vertex -> label
    max.summary <- final.summary[final.summary$type == "max", ]
    min.summary <- final.summary[final.summary$type == "min", ]

    max.vertex.to.label <- setNames(max.summary$label, max.summary$vertex)
    min.vertex.to.label <- setNames(min.summary$label, min.summary$vertex)

    ## Name the basin list elements by looking up each basin's vertex
    if (length(current.basins$lmax_basins) > 0) {
        basin.vertices <- sapply(current.basins$lmax_basins, function(b) b$vertex)
        names(current.basins$lmax_basins) <- as.character(max.vertex.to.label[as.character(basin.vertices)])
    }

    if (length(current.basins$lmin_basins) > 0) {
        basin.vertices <- sapply(current.basins$lmin_basins, function(b) b$vertex)
        names(current.basins$lmin_basins) <- as.character(min.vertex.to.label[as.character(basin.vertices)])
    }

    if (verbose) {
        cat("  Assigned", length(current.basins$lmax_basins), "maximum labels\n")
        cat("  Assigned", length(current.basins$lmin_basins), "minimum labels\n")
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

    if (verbose) {
        n.max.final <- sum(final.summary$type == "max")
        n.min.final <- sum(final.summary$type == "min")
        cat("\nRefinement complete!\n")
        cat(sprintf("Final structure: %d maxima and %d minima\n", n.max.final, n.min.final))
    }

    ## Compute final overlap distance matrices for diagnostics
    if (verbose) {
        cat("Computing final overlap distance matrices...\n")
    }

    max.overlap.dist <- NULL
    min.overlap.dist <- NULL

    ## Compute overlap distances for maxima
    if (length(current.basins$lmax_basins) > 1) {
        max.summary <- final.summary[final.summary$type == "max", ]
        max.vertices.list <- list()
        for (i in seq_len(nrow(max.summary))) {
            label <- max.summary$label[i]
            vertex <- max.summary$vertex[i]
            for (basin in current.basins$lmax_basins) {
                if (basin$vertex == vertex) {
                    max.vertices.list[[label]] <- basin$basin_df[, 1]
                    break
                }
            }
        }
        max.overlap.dist <- compute.overlap.distance.matrix(max.vertices.list)
    }

    ## Compute overlap distances for minima
    if (length(current.basins$lmin_basins) > 1) {
        min.summary <- final.summary[final.summary$type == "min", ]
        min.vertices.list <- list()
        for (i in seq_len(nrow(min.summary))) {
            label <- min.summary$label[i]
            vertex <- min.summary$vertex[i]
            for (basin in current.basins$lmin_basins) {
                if (basin$vertex == vertex) {
                    min.vertices.list[[label]] <- basin$basin_df[, 1]
                    break
                }
            }
        }
        min.overlap.dist <- compute.overlap.distance.matrix(min.vertices.list)
    }

    ## Return both basins and summary
    result <- list(
        basins = current.basins,
        summary = final.summary,
        max.overlap.dist = max.overlap.dist,
        min.overlap.dist = min.overlap.dist
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

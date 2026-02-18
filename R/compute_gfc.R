#' Compute Gradient Flow Complex
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
#'   containing fewer than this many vertices are removed. Default is 10.
#' @param expand.basins Logical indicating whether to expand basins to cover all
#'   graph vertices (default: TRUE). When TRUE, vertices not covered by any
#'   retained basin are assigned to the nearest basin based on shortest path
#'   distance in the weighted graph. This ensures complete coverage after
#'   filtering removes spurious extrema.
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
#'     \item{max.vertices.list}{Named list of integer vectors containing vertex
#'       indices for each maximum basin. Names match basin labels from the summary.}
#'     \item{min.vertices.list}{Named list of integer vectors containing vertex
#'       indices for each minimum basin. Names match basin labels from the summary.}
#'     \item{expanded.max.vertices.list}{If \code{expand.basins = TRUE}, a named
#'       list of expanded maximum basins covering all graph vertices. Otherwise
#'       \code{NULL}.}
#'     \item{expanded.min.vertices.list}{If \code{expand.basins = TRUE}, a named
#'       list of expanded minimum basins covering all graph vertices. Otherwise
#'       \code{NULL}.}
#'     \item{stage.history}{A data frame recording the number of maxima and minima
#'       before and after each refinement stage, facilitating reproducible reporting.}
#'     \item{parameters}{A named list of all parameter values used in the refinement
#'       process, enabling exact replication and automated report generation.}
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
#' determine which basins to retain based on these geometric measures. Additionally,
#' basins with fewer than \code{min.basin.size} vertices are removed as they lack
#' sufficient statistical support.
#'
#' When \code{expand.basins = TRUE}, an additional step assigns uncovered vertices
#' to their nearest retained basin. After filtering spurious extrema, some vertices
#' may not belong to any retained basin. These are reassigned based on shortest
#' path distance in the weighted graph, ensuring complete coverage for downstream
#' analyses that require every vertex to have a basin assignment.
#'
#' Each filtering stage can be disabled by setting the corresponding logical
#' parameter to FALSE, allowing users to customize the refinement pipeline
#' for specific applications.
#'
#' @examples
#' \dontrun{
#' # Compute refined basins with default parameters
#' result <- compute.gfc(adj.list,
#'                                  edge.length.list,
#'                                  fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  verbose = TRUE)
#'
#' # Access the refined basins and summary
#' refined.basins <- result$basins
#' basin.summary <- result$summary
#'
#' # View stage history
#' print(result$stage.history)
#'
#' # Generate methods paragraph
#' report <- generate.refinement.report(result)
#' cat(report)
#'
#' # Use custom thresholds for more aggressive filtering
#' result <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  min.rel.value.max = 1.2,
#'                                  max.rel.value.min = 0.8,
#'                                  min.basin.size = 10,
#'                                  expand.basins = TRUE,
#'                                  p.mean.hopk.dist.threshold = 0.85,
#'                                  verbose = TRUE)
#'
#' # Skip clustering steps
#' result <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  apply.maxima.clustering = FALSE,
#'                                  apply.minima.clustering = FALSE,
#'                                  verbose = TRUE)
#' }
#'
#' @seealso \code{\link{generate.refinement.report}} for automated report generation,
#'   \code{\link{expand.basins.to.cover}} for the basin expansion algorithm
#'
#' @export
compute.gfc <- function(adj.list,
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
                                   min.basin.size = 10,
                                   expand.basins = TRUE,
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

    ## Initialize stage history tracking
    stage.history <- data.frame(
        stage = character(),
        description = character(),
        n.max.before = integer(),
        n.max.after = integer(),
        n.min.before = integer(),
        n.min.after = integer(),
        stringsAsFactors = FALSE
    )

    ## -------------------------------------------------------
    ## Step 0: Breaking ties (if any)
    ## -------------------------------------------------------
    fitted.values <- break.ties(fitted.values,
                                noise.scale = 1e-10,
                                min.abs.noise = 1e-12,
                                preserve.bounds = TRUE,
                                seed = NULL,
                                verbose = FALSE)

    ## -------------------------------------------------------
    ## Step 1: Compute initial basins of attraction
    ## -------------------------------------------------------
    if (verbose) {
        cat("Step 1: Computing initial basins of attraction...\n")
    }

    current.basins <- compute.basins.of.attraction(adj.list,
                                                   edge.length.list,
                                                   fitted.values,
                                                   edge.length.quantile.thld,
                                                   with.trajectories)

    summary.use.cpp <- TRUE
    vertex.metrics <- NULL
    vertex.metrics <- tryCatch({
        if (verbose) {
            cat("Precomputing vertex metrics for basin summaries (C++)...\n")
        }
        adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1L))
        .Call(
            "S_precompute_basin_vertex_metrics",
            adj.list.0based,
            edge.length.list,
            hop.k,
            PACKAGE = "gflow"
        )
    }, error = function(e) {
        summary.use.cpp <<- FALSE
        warning(
            "C++ basin-summary metric precompute failed; using R summary fallback. ",
            "Error: ", conditionMessage(e),
            call. = FALSE
        )
        NULL
    })

    initial.summary <- summary(current.basins,
                               adj.list,
                               edge.length.list,
                               hop.k,
                               vertex.metrics = vertex.metrics,
                               use.cpp = summary.use.cpp)
    n.max <- sum(initial.summary$type == "max")
    n.min <- sum(initial.summary$type == "min")

    stage.history <- rbind(stage.history, data.frame(
        stage = "initial",
        description = sprintf("Initial basins (edge quantile = %.2f)",
                              edge.length.quantile.thld),
        n.max.before = NA_integer_,
        n.max.after = n.max,
        n.min.before = NA_integer_,
        n.min.after = n.min,
        stringsAsFactors = FALSE
    ))

    if (verbose) {
        cat(sprintf("  Found %d maxima and %d minima\n", n.max, n.min))
        cat("initial.summary:\n")
        print(initial.summary)
    }

    ## -------------------------------------------------------
    ## Step 2: Filter by relative values
    ## -------------------------------------------------------
    if (apply.relvalue.filter) {
        current.summary <- summary(current.basins,
                                   adj.list,
                                   edge.length.list,
                                   hop.k,
                                   vertex.metrics = vertex.metrics,
                                   use.cpp = summary.use.cpp)
        n.max.before <- sum(current.summary$type == "max")
        n.min.before <- sum(current.summary$type == "min")

        if (verbose) {
            cat(sprintf("\nStep 2: Filtering by relative values (max >= %.2f, min <= %.2f)...\n",
                        min.rel.value.max, max.rel.value.min))
        }

        current.basins <- filter.basins.by.relvalue(current.basins,
                                                    min.rel.value.max = min.rel.value.max,
                                                    max.rel.value.min = max.rel.value.min)

        relvalue.summary <- summary(current.basins,
                                    adj.list,
                                    edge.length.list,
                                    hop.k,
                                    vertex.metrics = vertex.metrics,
                                    use.cpp = summary.use.cpp)
        n.max.after <- sum(relvalue.summary$type == "max")
        n.min.after <- sum(relvalue.summary$type == "min")

        stage.history <- rbind(stage.history, data.frame(
            stage = "relvalue",
            description = sprintf("Relative value filter (max >= %.2f, min <= %.2f)",
                                  min.rel.value.max, max.rel.value.min),
            n.max.before = n.max.before,
            n.max.after = n.max.after,
            n.min.before = n.min.before,
            n.min.after = n.min.after,
            stringsAsFactors = FALSE
        ))

        if (verbose) {
            cat(sprintf("  Retained %d maxima and %d minima\n", n.max.after, n.min.after))
            cat("relvalue.summary:\n")
            print(relvalue.summary)
        }
    } else {
        if (verbose) {
            cat("\nStep 2: Skipping relative value filtering\n")
        }
    }

    ## -------------------------------------------------------
    ## Step 3: Cluster and merge maxima
    ## -------------------------------------------------------
    if (apply.maxima.clustering) {
        current.summary <- summary(current.basins,
                                   adj.list,
                                   edge.length.list,
                                   hop.k,
                                   vertex.metrics = vertex.metrics,
                                   use.cpp = summary.use.cpp)
        n.max.before <- sum(current.summary$type == "max")
        n.min.before <- sum(current.summary$type == "min")

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
            n.clusters <- length(unique(max.clusters$cluster.assignments))
            n.singletons <- sum(table(max.clusters$cluster.assignments) == 1)
            n.merged <- n.clusters - n.singletons
            cat(sprintf("  Found %d clusters (%d singletons, %d to be merged)\n",
                        n.clusters, n.singletons, n.merged))
            cat("  Merging clustered maxima...\n")
        }

        current.basins <- merge.clustered.basins(current.basins,
                                                 max.clusters,
                                                 extrema.type = "max")

        merged.max.summary <- summary(current.basins,
                                      adj.list,
                                      edge.length.list,
                                      hop.k,
                                      vertex.metrics = vertex.metrics,
                                      use.cpp = summary.use.cpp)
        n.max.after <- sum(merged.max.summary$type == "max")
        n.min.after <- sum(merged.max.summary$type == "min")

        stage.history <- rbind(stage.history, data.frame(
            stage = "merge.max",
            description = sprintf("Cluster and merge maxima (overlap threshold = %.2f)",
                                  max.overlap.threshold),
            n.max.before = n.max.before,
            n.max.after = n.max.after,
            n.min.before = n.min.before,
            n.min.after = n.min.after,
            stringsAsFactors = FALSE
        ))

        if (verbose) {
            cat(sprintf("  Result: %d maxima after merging\n", n.max.after))
            cat("merged.max.summary:\n")
            print(merged.max.summary)
        }
    } else {
        if (verbose) {
            cat("\nStep 3: Skipping maxima clustering\n")
        }
    }

    ## -------------------------------------------------------
    ## Step 4: Cluster and merge minima
    ## -------------------------------------------------------
    if (apply.minima.clustering) {
        current.summary <- summary(current.basins,
                                   adj.list,
                                   edge.length.list,
                                   hop.k,
                                   vertex.metrics = vertex.metrics,
                                   use.cpp = summary.use.cpp)
        n.max.before <- sum(current.summary$type == "max")
        n.min.before <- sum(current.summary$type == "min")

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
            n.clusters <- length(unique(min.clusters$cluster.assignments))
            n.singletons <- sum(table(min.clusters$cluster.assignments) == 1)
            n.merged <- n.clusters - n.singletons
            cat(sprintf("  Found %d clusters (%d singletons, %d to be merged)\n",
                        n.clusters, n.singletons, n.merged))
            cat("  Merging clustered minima...\n")
        }

        current.basins <- merge.clustered.basins(current.basins,
                                                 min.clusters,
                                                 extrema.type = "min")

        merged.min.summary <- summary(current.basins,
                                      adj.list,
                                      edge.length.list,
                                      hop.k,
                                      vertex.metrics = vertex.metrics,
                                      use.cpp = summary.use.cpp)
        n.max.after <- sum(merged.min.summary$type == "max")
        n.min.after <- sum(merged.min.summary$type == "min")

        stage.history <- rbind(stage.history, data.frame(
            stage = "merge.min",
            description = sprintf("Cluster and merge minima (overlap threshold = %.2f)",
                                  min.overlap.threshold),
            n.max.before = n.max.before,
            n.max.after = n.max.after,
            n.min.before = n.min.before,
            n.min.after = n.min.after,
            stringsAsFactors = FALSE
        ))

        if (verbose) {
            cat(sprintf("  Result: %d minima after merging\n", n.min.after))
            cat("merged.min.summary:\n")
            print(merged.min.summary)
        }
    } else {
        if (verbose) {
            cat("\nStep 4: Skipping minima clustering\n")
        }
    }

    ## ------------------------------------------------------------
    ## Step 5: Filter by geometric characteristics and basin size
    ## ------------------------------------------------------------
    if (apply.geometric.filter) {
        current.summary <- summary(current.basins,
                                   adj.list,
                                   edge.length.list,
                                   hop.k,
                                   vertex.metrics = vertex.metrics,
                                   use.cpp = summary.use.cpp)
        n.max.before <- sum(current.summary$type == "max")
        n.min.before <- sum(current.summary$type == "min")

        if (verbose) {
            cat(sprintf("\nStep 5: Filtering by geometric characteristics and basin size:\n"))
            cat(sprintf("  p.mean.nbrs.dist < %.2f, p.mean.hopk.dist < %.2f, basin.size >= %d ...\n",
                        p.mean.nbrs.dist.threshold, p.mean.hopk.dist.threshold,
                        min.basin.size))
        }

        ## Get summary with geometric characteristics
        geom.summary <- summary(current.basins,
                                adj.list,
                                edge.length.list,
                                hop.k = hop.k,
                                vertex.metrics = vertex.metrics,
                                use.cpp = summary.use.cpp)

        ## Filter maxima
        max.basin.df <- geom.summary[geom.summary$type == "max", ]
        good.max.vertices <- max.basin.df$vertex[
                                              max.basin.df$p.mean.nbrs.dist < p.mean.nbrs.dist.threshold &
                                              max.basin.df$p.mean.hopk.dist < p.mean.hopk.dist.threshold &
                                              max.basin.df$basin.size >= min.basin.size
                                          ]

        if (verbose) {
            n.max.after <- length(good.max.vertices)
            cat(sprintf("  Maxima: %d -> %d (removed %d)\n",
                        n.max.before, n.max.after, n.max.before - n.max.after))
            cat("geom.summary:\n")
            print(geom.summary)
        }

        ## Filter minima
        min.basin.df <- geom.summary[geom.summary$type == "min", ]
        good.min.vertices <- min.basin.df$vertex[
                                              min.basin.df$p.mean.nbrs.dist < p.mean.nbrs.dist.threshold &
                                              min.basin.df$p.mean.hopk.dist < p.mean.hopk.dist.threshold &
                                              min.basin.df$basin.size >= min.basin.size
                                          ]

        if (verbose) {
            n.min.after <- length(good.min.vertices)
            cat(sprintf("  Minima: %d -> %d (removed %d)\n",
                        n.min.before, n.min.after, n.min.before - n.min.after))
        }

        ## Combine good extrema and filter basins
        good.extrema <- c(good.min.vertices, good.max.vertices)
        if (length(good.extrema) > 0) {
            current.basins <- filter.basins(current.basins, good.extrema)
        } else {
            if (verbose) {
                cat("  No extrema passed geometric filter; retaining empty basin set\n")
            }
            current.basins$lmin_basins <- list()
            current.basins$lmax_basins <- list()
        }

        n.max.after <- length(good.max.vertices)
        n.min.after <- length(good.min.vertices)

        stage.history <- rbind(stage.history, data.frame(
            stage = "geometric",
            description = sprintf("Geometric filter (p.mean.nbrs.dist < %.2f, p.mean.hopk.dist < %.2f, size >= %d)",
                                  p.mean.nbrs.dist.threshold, p.mean.hopk.dist.threshold, min.basin.size),
            n.max.before = n.max.before,
            n.max.after = n.max.after,
            n.min.before = n.min.before,
            n.min.after = n.min.after,
            stringsAsFactors = FALSE
        ))

    } else {
        if (verbose) {
            cat("\nStep 5: Skipping geometric filtering\n")
        }
    }

    ## -------------------------------------------------------
    ## Generate final summary
    ## -------------------------------------------------------
    if (verbose) {
        cat(sprintf("\nGenerating final summary (hop.k = %d)...\n", hop.k))
    }

    final.summary <- summary(current.basins,
                             adj.list,
                             edge.length.list,
                             hop.k = hop.k,
                             vertex.metrics = vertex.metrics,
                             use.cpp = summary.use.cpp)

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
    if (verbose) {
        cat("Storing graph structure in basins object...\n")
    }

    if (is.null(current.basins$adj.list)) {
        current.basins$adj.list <- adj.list
    }

    if (is.null(current.basins$edge.length.list)) {
        current.basins$edge.length.list <- edge.length.list
    }

    current.basins$hop.k <- hop.k

    ## Construct vertices lists for maxima and minima
    if (verbose) {
        cat("Constructing basin vertices lists...\n")
    }

    max.vertices.list <- list()
    if (nrow(max.summary) > 0) {
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
    }

    min.vertices.list <- list()
    if (nrow(min.summary) > 0) {
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
    }

    ## ----------------------------------------------------------
    ## Compute final overlap distance matrices for diagnostics
    ## ----------------------------------------------------------
    if (verbose) {
        cat("Computing final overlap distance matrices...\n")
    }

    max.overlap.dist <- NULL
    min.overlap.dist <- NULL

    if (length(max.vertices.list) > 1) {
        max.overlap.dist <- compute.overlap.distance.matrix(max.vertices.list)
    }

    if (length(min.vertices.list) > 1) {
        min.overlap.dist <- compute.overlap.distance.matrix(min.vertices.list)
    }

    ## Step 6 (optional): Expand basins to cover all vertices
    expanded.max.vertices.list <- NULL
    expanded.min.vertices.list <- NULL

    if (expand.basins) {
        if (verbose) {
            cat("\nStep 6: Expanding basins to cover all vertices...\n")
        }

        n.vertices <- length(adj.list)

        ## Count coverage before expansion
        n.covered.max <- length(unique(unlist(max.vertices.list)))
        n.covered.min <- length(unique(unlist(min.vertices.list)))

        ## Expand maximum basins
        if (length(max.vertices.list) > 0) {
            expanded.max.vertices.list <- expand.basins.to.cover(
                basins.vertices.list = max.vertices.list,
                adj.list = adj.list,
                weight.list = edge.length.list,
                n.vertices = n.vertices
            )
            n.expanded.max <- length(unique(unlist(expanded.max.vertices.list)))

            if (verbose) {
                cat(sprintf("  Maxima basins: %d -> %d vertices covered\n",
                            n.covered.max, n.expanded.max))
            }
        }

        ## Expand minimum basins
        if (length(min.vertices.list) > 0) {
            expanded.min.vertices.list <- expand.basins.to.cover(
                basins.vertices.list = min.vertices.list,
                adj.list = adj.list,
                weight.list = edge.length.list,
                n.vertices = n.vertices
            )
            n.expanded.min <- length(unique(unlist(expanded.min.vertices.list)))

            if (verbose) {
                cat(sprintf("  Minima basins: %d -> %d vertices covered\n",
                            n.covered.min, n.expanded.min))
            }
        }

        stage.history <- rbind(stage.history, data.frame(
            stage = "expand",
            description = "Expand basins to cover all vertices",
            n.max.before = n.covered.max,
            n.max.after = n.vertices,
            n.min.before = n.covered.min,
            n.min.after = n.vertices,
            stringsAsFactors = FALSE
        ))

    } else {
        if (verbose) {
            cat("\nStep 6: Skipping basin expansion\n")
        }
    }

    ## Store parameters for reporting
    parameters <- list(
        edge.length.quantile.thld = edge.length.quantile.thld,
        min.rel.value.max = min.rel.value.max,
        max.rel.value.min = max.rel.value.min,
        max.overlap.threshold = max.overlap.threshold,
        min.overlap.threshold = min.overlap.threshold,
        p.mean.nbrs.dist.threshold = p.mean.nbrs.dist.threshold,
        p.mean.hopk.dist.threshold = p.mean.hopk.dist.threshold,
        p.deg.threshold = p.deg.threshold,
        min.basin.size = min.basin.size,
        expand.basins = expand.basins,
        hop.k = hop.k
    )

    if (verbose) {
        n.max.final <- sum(final.summary$type == "max")
        n.min.final <- sum(final.summary$type == "min")
        cat("\nRefinement complete!\n")
        cat(sprintf("Final structure: %d maxima and %d minima\n", n.max.final, n.min.final))
    }

    ## Return basins, summary, diagnostics, history, and parameters
    result <- list(
        basins = current.basins,
        summary = final.summary,
        max.overlap.dist = max.overlap.dist,
        min.overlap.dist = min.overlap.dist,
        max.vertices.list = max.vertices.list,
        min.vertices.list = min.vertices.list,
        expanded.max.vertices.list = expanded.max.vertices.list,
        expanded.min.vertices.list = expanded.min.vertices.list,
        stage.history = stage.history,
        parameters = parameters
    )

    class(result) <- "basins_of_attraction"

    return(result)
}


#' Expand Basins to Cover All Graph Vertices
#'
#' Assigns uncovered vertices to their nearest basin based on shortest path
#' distance in the graph. This is useful after filtering spurious extrema,
#' when some vertices may no longer belong to any retained basin.
#'
#' @param basins.vertices.list A named list of integer vectors, where each
#'   element contains the vertex indices belonging to a basin. Names should
#'   correspond to basin labels (e.g., "M1", "M2", ... or "m1", "m2", ...).
#' @param adj.list Adjacency list representation of the graph. Element \code{i}
#'   contains the indices of vertices adjacent to vertex \code{i}.
#' @param weight.list List of edge weights (typically edge lengths) corresponding
#'   to the adjacency list. Element \code{i} contains the weights of edges
#'   incident to vertex \code{i}.
#' @param n.vertices Integer specifying the total number of vertices in the graph.
#'   If \code{NULL} (default), inferred from the length of \code{adj.list}.
#'
#' @return A named list with the same structure as \code{basins.vertices.list},
#'   where each basin has been expanded to include its nearest uncovered vertices.
#'   The union of all returned basins covers all graph vertices.
#'
#' @details
#' After filtering spurious local extrema, the union of retained basins may not
#' cover the entire graph. Vertices that belonged to removed basins are left
#' unassigned. This function reassigns such vertices to the nearest retained
#' basin based on shortest path distance in the weighted graph.
#'
#' For each uncovered vertex, the function computes the shortest path distance
#' to all covered vertices, identifies the closest covered vertex, and assigns
#' the uncovered vertex to the basin containing that closest vertex. Ties are
#' broken arbitrarily (by the order in which basins appear in the input list).
#'
#' The function uses \code{igraph::distances()} for efficient shortest path
#' computation, restricting the calculation to distances from uncovered vertices
#' to covered vertices rather than computing the full distance matrix.
#'
#' @examples
#' \dontrun{
#' # After computing refined basins
#' result <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  min.basin.size = 10)
#'
#' # Expand maximum basins to cover all vertices
#' expanded.max.basins <- expand.basins.to.cover(
#'     basins.vertices.list = result$max.vertices.list,
#'     adj.list = adj.list,
#'     weight.list = edge.length.list,
#'     n.vertices = length(adj.list)
#' )
#'
#' # Verify full coverage
#' all.vertices <- unique(unlist(expanded.max.basins))
#' length(all.vertices) == length(adj.list)  # Should be TRUE
#' }
#'
#' @seealso \code{\link{compute.gfc}} for computing refined basins
#'
#' @importFrom igraph graph_from_data_frame distances E
#'
#' @export
expand.basins.to.cover <- function(basins.vertices.list,
                                   adj.list,
                                   weight.list,
                                   n.vertices = NULL) {

    ## Validate inputs
    if (!is.list(basins.vertices.list) || length(basins.vertices.list) == 0) {
        stop("basins.vertices.list must be a non-empty list")
    }

    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("adj.list and weight.list must be lists")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (is.null(n.vertices)) {
        n.vertices <- length(adj.list)
    }

    ## Compute union of all basin vertices
    covered.vertices <- unique(unlist(basins.vertices.list))

    ## Find uncovered vertices
    all.vertices <- seq_len(n.vertices)
    uncovered.vertices <- setdiff(all.vertices, covered.vertices)

    ## If all vertices are covered, return input unchanged
    if (length(uncovered.vertices) == 0) {
        return(basins.vertices.list)
    }

    ## Create lookup: covered vertex -> basin label
    vertex.to.basin <- character(n.vertices)
    for (basin.label in names(basins.vertices.list)) {
        basin.vertices <- basins.vertices.list[[basin.label]]
        vertex.to.basin[basin.vertices] <- basin.label
    }

    ## Build igraph from adjacency list and weight list
    edge.list <- data.frame(
        from = integer(0),
        to = integer(0),
        weight = numeric(0)
    )

    for (i in seq_len(n.vertices)) {
        if (length(adj.list[[i]]) > 0) {
            neighbors <- adj.list[[i]]
            ## Only add edges in one direction (i < j) to avoid duplicates
            mask <- neighbors > i
            if (any(mask)) {
                new.edges <- data.frame(
                    from = i,
                    to = neighbors[mask],
                    weight = weight.list[[i]][mask]
                )
                edge.list <- rbind(edge.list, new.edges)
            }
        }
    }

    ## Create the graph
    g <- igraph::graph_from_data_frame(
        edge.list,
        directed = FALSE,
        vertices = data.frame(name = as.character(seq_len(n.vertices)))
    )

    ## Compute distances from uncovered vertices to covered vertices only
    dist.matrix <- igraph::distances(
        g,
        v = as.character(uncovered.vertices),
        to = as.character(covered.vertices),
        weights = igraph::E(g)$weight
    )

    ## Initialize expanded basins as copies of input
    expanded.basins <- lapply(basins.vertices.list, function(x) x)

    ## Assign each uncovered vertex to its nearest basin
    for (i in seq_along(uncovered.vertices)) {
        v <- uncovered.vertices[i]
        distances.to.covered <- dist.matrix[i, ]

        ## Find closest covered vertex
        min.idx <- which.min(distances.to.covered)
        closest.covered <- covered.vertices[min.idx]

        ## Look up which basin the closest covered vertex belongs to
        target.basin <- vertex.to.basin[closest.covered]

        ## Add uncovered vertex to that basin
        expanded.basins[[target.basin]] <- c(expanded.basins[[target.basin]], v)
    }

    return(expanded.basins)
}


#' Generate Text Report of Basin Refinement Process
#'
#' Creates a paragraph suitable for inclusion in a methods section describing
#' the basin refinement process with specific parameter values and counts.
#'
#' @param basins.result Object returned by \code{compute.gfc()}
#' @param response.name Character string naming the response variable in LaTeX
#'   format, for example \code{"\\widehat{sPTB}"}.
#'
#' @return Character string containing the formatted paragraph ready for
#'   inclusion in a LaTeX document.
#'
#' @examples
#' \dontrun{
#' result <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                                  edge.length.quantile.thld = 0.9,
#'                                  min.basin.size = 10,
#'                                  expand.basins = TRUE,
#'                                  verbose = TRUE)
#'
#' report <- generate.refinement.report(result)
#' cat(report)
#'
#' ## Write to file for inclusion in Rmd/LaTeX document
#' writeLines(report, "refinement_methods.tex")
#' }
#'
#' @seealso \code{\link{compute.gfc}} for computing refined basins
#'
#' @export
generate.refinement.report <- function(basins.result,
                                       response.name = "\\widehat{sPTB}") {

    params <- basins.result$parameters
    history <- basins.result$stage.history

    ## Extract counts from each stage
    initial <- history[history$stage == "initial", ]
    relvalue <- history[history$stage == "relvalue", ]
    merge.max <- history[history$stage == "merge.max", ]
    merge.min <- history[history$stage == "merge.min", ]
    geometric <- history[history$stage == "geometric", ]
    expand <- history[history$stage == "expand", ]

    ## Handle case where some stages may have been skipped
    if (nrow(relvalue) == 0) {
        relvalue <- initial
        relvalue$n.max.before <- initial$n.max.after
        relvalue$n.min.before <- initial$n.min.after
    }

    if (nrow(merge.max) == 0) {
        merge.max <- relvalue
        merge.max$n.max.before <- relvalue$n.max.after
        merge.max$n.min.before <- relvalue$n.min.after
    }

    if (nrow(merge.min) == 0) {
        merge.min <- merge.max
        merge.min$n.max.before <- merge.max$n.max.after
        merge.min$n.min.before <- merge.max$n.min.after
    }

    if (nrow(geometric) == 0) {
        geometric <- merge.min
        geometric$n.max.before <- merge.min$n.max.after
        geometric$n.min.before <- merge.min$n.min.after
    }

    ## Build the main paragraph
    para <- sprintf(
"The refinement begins by computing basins of attraction for all local extrema of $%s$ on the ikNN graph. To prevent trajectories from jumping across distant regions of the graph, we restricted basin computation to edges below the %.0fth percentile of the edge length distribution. This initial step identified %d local maxima and %d local minima.

The first filtering stage removes extrema whose values lie too close to the median of $%s$. We define the relative value of an extremum as its fitted value divided by the median. Maxima with relative values below %.2f and minima with relative values above %.2f were eliminated, as such extrema represent minor fluctuations rather than prominent features of the response surface. This step retained %d maxima and %d minima.

The second stage addresses geometric redundancy by clustering extrema whose basins exhibit substantial overlap. We applied an overlap threshold of %.0f%% for maxima and %.0f%% for minima. Within each cluster, we merged the constituent basins, retaining the extremum with the most extreme value as the representative. This clustering and merging process reduced the count to %d maxima and %d minima.

The final automated stage filters extrema based on structural characteristics of their basins and basin size. We removed extrema exceeding the %.0fth percentile threshold for mean distance to neighbors or degree, and eliminated extrema whose basins contained fewer than %d samples. This filtering yielded %d maxima and %d minima.",
        response.name,
        params$edge.length.quantile.thld * 100,
        initial$n.max.after, initial$n.min.after,
        response.name,
        params$min.rel.value.max, params$max.rel.value.min,
        relvalue$n.max.after, relvalue$n.min.after,
        params$max.overlap.threshold * 100, params$min.overlap.threshold * 100,
        merge.min$n.max.after, merge.min$n.min.after,
        params$p.mean.hopk.dist.threshold * 100, params$min.basin.size,
        geometric$n.max.after, geometric$n.min.after
    )

    ## Add expansion paragraph if applicable
    if (nrow(expand) > 0 && params$expand.basins) {
        n.vertices <- expand$n.max.after  # Total vertices after expansion
        n.uncovered.max <- n.vertices - expand$n.max.before
        n.uncovered.min <- n.vertices - expand$n.min.before

        expand.para <- sprintf(
"

After filtering, %d vertices were not covered by any maximum basin and %d vertices were not covered by any minimum basin. These uncovered vertices, which originally belonged to basins of removed spurious extrema, were reassigned to the nearest retained basin based on shortest path distance in the weighted graph. This ensures complete coverage of all %d graph vertices for downstream analyses.",
            n.uncovered.max, n.uncovered.min, n.vertices
        )

        para <- paste0(para, expand.para)
    }

    return(para)
}


#' Print Method for Refined Basins
#'
#' @param x A basins_of_attraction object from compute.gfc()
#' @param ... Additional arguments passed to print methods
#'
#' @export
print.basins_of_attraction <- function(x, ...) {
    cat("Refined Basins of Attraction\n")
    cat("============================\n\n")

    n.max <- sum(x$summary$type == "max")
    n.min <- sum(x$summary$type == "min")

    cat(sprintf("Number of maxima: %d\n", n.max))
    cat(sprintf("Number of minima: %d\n", n.min))
    cat(sprintf("Total extrema: %d\n\n", n.max + n.min))

    ## Report coverage
    if (!is.null(x$max.vertices.list)) {
        n.covered.max <- length(unique(unlist(x$max.vertices.list)))
        cat(sprintf("Vertices covered by maxima basins: %d\n", n.covered.max))
    }
    if (!is.null(x$min.vertices.list)) {
        n.covered.min <- length(unique(unlist(x$min.vertices.list)))
        cat(sprintf("Vertices covered by minima basins: %d\n", n.covered.min))
    }

    ## Report expansion if applicable
    if (!is.null(x$expanded.max.vertices.list)) {
        n.expanded <- length(unique(unlist(x$expanded.max.vertices.list)))
        cat(sprintf("Vertices after maxima expansion: %d\n", n.expanded))
    }
    if (!is.null(x$expanded.min.vertices.list)) {
        n.expanded <- length(unique(unlist(x$expanded.min.vertices.list)))
        cat(sprintf("Vertices after minima expansion: %d\n", n.expanded))
    }

    cat("\nBasin summary (first 10 rows):\n")
    print(head(x$summary, 10))

    if (nrow(x$summary) > 10) {
        cat(sprintf("\n... and %d more rows\n", nrow(x$summary) - 10))
    }

    invisible(x)
}

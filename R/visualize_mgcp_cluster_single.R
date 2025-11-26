#' Visualize a Single MGCP Cluster
#'
#' Creates a detailed heatmap visualization of a single cluster from MGCP results,
#' showing all vertices in that cluster across all phylotypes (features).
#'
#' @param data.matrix Numeric matrix to be visualized (n.vertices Ã— n.phylotypes).
#' @param mgcp.result Result object from \code{extract.stable.clusters()}.
#' @param cluster.id Integer specifying which cluster to visualize (1, 2, 3, ...).
#' @param cluster.rows Logical indicating whether to perform hierarchical clustering
#'   on rows within the cluster. Default is \code{TRUE} to reveal internal structure.
#' @param cluster.cols Logical indicating whether to perform hierarchical clustering
#'   on columns. Default is \code{TRUE}.
#' @param show.phylotype.names Logical indicating whether to display phylotype
#'   names. Default is \code{TRUE}.
#' @param show.statistics Logical indicating whether to display cluster statistics
#'   (size, gap, rank structure) as text annotation. Default is \code{TRUE}.
#' @param main.title Character string for heatmap title. If \code{NULL}, generates
#'   automatic title with cluster information.
#' @param color.scheme Character specifying color scheme: "default" (blue-white-red),
#'   "viridis", "plasma", or "red-blue". Default is "default".
#' @param phylotype.fontsize Numeric specifying font size for phylotype names.
#'   Default is 8.
#' @param verbose Logical indicating whether to print summary statistics and
#'   interpretation to console. Default is \code{TRUE}.
#' @param phylo.mean.thld1 Numeric threshold for moderate phylotype signal strength.
#'   Phylotypes with absolute mean co-monotonicity between \code{phylo.mean.thld1}
#'   and \code{phylo.mean.thld2} are classified as "moderate". Default is 0.1.
#' @param phylo.mean.thld2 Numeric threshold for strong phylotype signal strength.
#'   Phylotypes with absolute mean co-monotonicity above \code{phylo.mean.thld2}
#'   are classified as "strong". Default is 0.3.
#' @param phylo.filter.thld Numeric threshold for phylotype filtering. Only phylotypes
#'   with absolute mean co-monotonicity of at least \code{phylo.filter.thld} are
#'   retained in visualization. Must be >= 0 and < 1. Default is 0 (no filtering).
#'   Higher values focus visualization on phylotypes with stronger signals.
#'
#' @return Returns invisibly a list containing heatmap object, cluster data, and
#'   statistics.
#'
#' @details
#' This function creates a detailed view of a single MGCP cluster, showing:
#' \itemize{
#'   \item All vertices assigned to the specified cluster
#'   \item All phylotypes (no filtering - MGCP philosophy)
#'   \item Cluster statistics (gap, size, rank structure)
#'   \item Optional hierarchical clustering to reveal internal patterns
#' }
#'
#' Unlike bi-clustering approaches that select feature subsets, this function
#' always displays all phylotypes, allowing you to see the complete co-monotonicity
#' pattern for the cluster without cherry-picking.
#'
#' @seealso \code{\link{visualize.mgcp.clusters}}, \code{\link{extract.stable.clusters}},
#'   \code{\link{analyze.extracted.clusters}}
#'
#' @examples
#' \dontrun{
#' # Extract clusters
#' mgcp.result <- extract.stable.clusters(comono.smooth, min.cluster.size = 20)
#'
#' # Visualize most stable cluster
#' viz <- visualize.mgcp.cluster(comono.smooth, mgcp.result, cluster.id = 1)
#'
#' # Visualize without clustering (preserve extraction order)
#' viz <- visualize.mgcp.cluster(
#'   comono.smooth,
#'   mgcp.result,
#'   cluster.id = 1,
#'   cluster.rows = FALSE
#' )
#'
#' # Hide phylotype names for cleaner view
#' viz <- visualize.mgcp.cluster(
#'   comono.smooth,
#'   mgcp.result,
#'   cluster.id = 4,
#'   show.phylotype.names = FALSE
#' )
#' }
#'
#' @export
visualize.mgcp.cluster <- function(data.matrix,
                                   mgcp.result,
                                   cluster.id,
                                   cluster.rows = TRUE,
                                   cluster.cols = TRUE,
                                   show.phylotype.names = TRUE,
                                   show.statistics = TRUE,
                                   main.title = NULL,
                                   color.scheme = c("default", "viridis", "plasma", "red-blue"),
                                   phylotype.fontsize = 8,
                                   verbose = TRUE,
                                   phylo.mean.thld1 = 0.1,
                                   phylo.mean.thld2 = 0.3,
                                   phylo.filter.thld = 0) {

  ## Validate phylotype threshold parameters
  if (!is.numeric(phylo.mean.thld1) || phylo.mean.thld1 < 0 || phylo.mean.thld1 >= 1) {
    stop("phylo.mean.thld1 must be a numeric value between 0 and 1")
  }
  
  if (!is.numeric(phylo.mean.thld2) || phylo.mean.thld2 < 0 || phylo.mean.thld2 >= 1) {
    stop("phylo.mean.thld2 must be a numeric value between 0 and 1")
  }
  
  if (phylo.mean.thld1 >= phylo.mean.thld2) {
    stop("phylo.mean.thld1 must be less than phylo.mean.thld2")
  }
  
  if (!is.numeric(phylo.filter.thld) || phylo.filter.thld < 0 || phylo.filter.thld >= 1) {
    stop("phylo.filter.thld must be a numeric value between 0 and 1 (exclusive)")
  }

  color.scheme <- match.arg(color.scheme)

    ## Validate inputs
    if (!is.numeric(cluster.id) || length(cluster.id) != 1) {
        stop("cluster.id must be a single integer")
    }

    cluster.id <- as.integer(cluster.id)

    if (cluster.id < 1 || cluster.id > mgcp.result$n.clusters) {
        stop(sprintf("cluster.id must be between 1 and %d (number of clusters found)",
                     mgcp.result$n.clusters))
    }

    ## Extract cluster assignments
    cluster.assignments <- mgcp.result$clusters
    cluster.members <- which(cluster.assignments == cluster.id)

    if (length(cluster.members) == 0) {
        stop(sprintf("No vertices found in cluster %d", cluster.id))
    }

    ## Extract cluster data
    cluster.data <- data.matrix[cluster.members, , drop = FALSE]
    
    ## Apply phylotype filtering if requested
    n.phylotypes.original <- ncol(cluster.data)
    if (phylo.filter.thld > 0) {
        phylo.means.prefilter <- colMeans(cluster.data)
        phylo.keep <- abs(phylo.means.prefilter) >= phylo.filter.thld
        
        if (sum(phylo.keep) == 0) {
            stop(sprintf("No phylotypes pass filter threshold %.3f. Lower phylo.filter.thld or choose different cluster.",
                        phylo.filter.thld))
        }
        
        cluster.data <- cluster.data[, phylo.keep, drop = FALSE]
        
        if (verbose) {
            cat(sprintf("Phylotype filtering (threshold = %.3f):\n", phylo.filter.thld))
            cat("  Original phylotypes:", n.phylotypes.original, "\n")
            cat("  Retained phylotypes:", ncol(cluster.data), "\n")
            cat("  Filtered out:", n.phylotypes.original - ncol(cluster.data), "\n\n")
        }
    }

    ## Get cluster statistics
    cluster.info <- mgcp.result$summary.stats[mgcp.result$summary.stats$cluster.id == cluster.id, ]

    if (nrow(cluster.info) == 0) {
        stop(sprintf("No statistics found for cluster %d", cluster.id))
    }

    ## Get detailed analysis if available
    cluster.analysis <- tryCatch({
        analyze.extracted.clusters(data.matrix, mgcp.result)
    }, error = function(e) {
        NULL
    })

    if (!is.null(cluster.analysis)) {
        cluster.stats <- cluster.analysis[cluster.analysis$cluster.id == cluster.id, ]
    } else {
        cluster.stats <- NULL
    }

    ## Print summary (if verbose)
    if (verbose) {
        cat("\n")
        cat("Cluster", cluster.id, "Summary:\n")
        cat(strrep("=", 50), "\n")
        cat("Size:", cluster.info$cluster.size, "vertices\n")
        cat("Gap (stability):", round(cluster.info$gap.size, 2), "\n")
        cat("Height:", round(cluster.info$cluster.height, 3), "\n")

        if (!is.null(cluster.stats)) {
            cat("Rank-1 variance:", round(cluster.stats$rank1.variance * 100, 1), "%\n")
            cat("Effective rank:", cluster.stats$effective.rank, "\n")
            cat("Median correlation:", round(cluster.stats$median.correlation, 3), "\n")
        }

        cat("Phylotypes displayed:", ncol(cluster.data), 
            if (phylo.filter.thld > 0) {
                sprintf("(filtered from %d)", n.phylotypes.original)
            } else {
                "(all retained)"
            }, "\n")
        cat("\n")
    }

    ## Identify characteristic phylotypes
    phylo.means <- colMeans(cluster.data)
    phylo.sds <- apply(cluster.data, 2, sd)

    strong.positive <- sum(phylo.means > phylo.mean.thld2)
    strong.negative <- sum(phylo.means < -phylo.mean.thld2)
    moderate <- sum(abs(phylo.means) >= phylo.mean.thld1 & abs(phylo.means) <= phylo.mean.thld2)
    weak <- sum(abs(phylo.means) < phylo.mean.thld1)

    if (verbose) {
        cat("Phylotype patterns:\n")
        cat(sprintf("  Strong positive (mean > %.2f): %d\n", phylo.mean.thld2, strong.positive))
        cat(sprintf("  Strong negative (mean < -%.2f): %d\n", phylo.mean.thld2, strong.negative))
        cat(sprintf("  Moderate (|mean| %.2f-%.2f): %d\n", phylo.mean.thld1, phylo.mean.thld2, moderate))
        cat(sprintf("  Weak (|mean| < %.2f): %d\n", phylo.mean.thld1, weak))
        cat("\n")
    }

    ## Show top phylotypes (if verbose)
    if (verbose) {
        if (strong.positive > 0) {
            top.positive <- head(sort(phylo.means[phylo.means > phylo.mean.thld2], decreasing = TRUE), 5)
            cat("Top positive phylotypes:\n")
            for (i in seq_along(top.positive)) {
                cat(sprintf("  %d. %s: %.3f\n", i, names(top.positive)[i], top.positive[i]))
            }
            cat("\n")
        }

        if (strong.negative > 0) {
            top.negative <- head(sort(phylo.means[phylo.means < -phylo.mean.thld2]), 5)
            cat("Top negative phylotypes:\n")
            for (i in seq_along(top.negative)) {
                cat(sprintf("  %d. %s: %.3f\n", i, names(top.negative)[i], top.negative[i]))
            }
            cat("\n")
        }
    }

    ## Set up color scheme
    if (color.scheme == "default") {
        col.fun <- circlize::colorRamp2(
            c(-1, -0.5, 0, 0.5, 1),
            c("blue", "lightblue", "white", "pink", "red")
        )
    } else if (color.scheme == "viridis") {
        col.fun <- circlize::colorRamp2(
            seq(-1, 1, length.out = 11),
            viridis::viridis(11)
        )
    } else if (color.scheme == "plasma") {
        col.fun <- circlize::colorRamp2(
            seq(-1, 1, length.out = 11),
            viridis::plasma(11)
        )
    } else if (color.scheme == "red-blue") {
        col.fun <- circlize::colorRamp2(
            c(-1, 0, 1),
            c("blue", "white", "red")
        )
    }

    ## Generate title with optional statistics
    if (is.null(main.title)) {
        if (show.statistics && !is.null(cluster.stats)) {

            main.title <- sprintf(
                "Cluster %d: %d vertices | Eff. rank: %d | Median cor: %.3f",
                cluster.id,
                cluster.info$cluster.size,
                cluster.stats$effective.rank,
                cluster.stats$median.correlation
            )

            ## main.title <- sprintf(
            ##     "Cluster %d: %d vertices (gap = %.1f) | Rank-1: %.1f%% | Eff. rank: %d | Median cor: %.3f",
            ##     cluster.id,
            ##     cluster.info$cluster.size,
            ##     cluster.info$gap.size,
            ##     cluster.stats$rank1.variance * 100,
            ##     cluster.stats$effective.rank,
            ##     cluster.stats$median.correlation
            ## )
        } else {
            main.title <- sprintf(
                "Cluster %d: %d vertices (gap = %.1f)",
                cluster.id,
                cluster.info$cluster.size,
                cluster.info$gap.size
            )
        }
    }

    ## No column annotation needed (statistics are in title)
    col.ha <- NULL

    ## Create row annotation showing phylotype means
    phylo.mean.colors <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("blue", "white", "red")
    )

    row.ha <- ComplexHeatmap::rowAnnotation(
        `Mean\nValue` = ComplexHeatmap::anno_barplot(
            phylo.means,
            gp = grid::gpar(fill = ifelse(phylo.means > 0, "red", "blue")),
            baseline = 0,
            axis = TRUE
        ),
        annotation_name_rot = 0,
        annotation_name_gp = grid::gpar(fontsize = 9),
        width = grid::unit(2, "cm")
    )

    ## Create heatmap
    ht <- ComplexHeatmap::Heatmap(
        t(cluster.data),  # Transpose: phylotypes as rows, vertices as columns
        name = "Co-mono",
        col = col.fun,
        cluster_rows = cluster.cols,  # Phylotypes
        cluster_columns = cluster.rows,  # Vertices
        show_row_names = show.phylotype.names,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = phylotype.fontsize),
        column_title = sprintf("Vertices (n = %d)", nrow(cluster.data)),
        column_title_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
        row_title = "Phylotypes",
        row_title_side = "left",
        row_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
        top_annotation = col.ha,
        right_annotation = row.ha,
        border = TRUE,
        heatmap_legend_param = list(
            title = "Co-mono",
            direction = "vertical",
            legend_height = grid::unit(4, "cm")
        ),
        use_raster = nrow(cluster.data) > 500 || ncol(cluster.data) > 500,
        raster_quality = 2
    )

    ## Draw heatmap
    ht.drawn <- ComplexHeatmap::draw(ht,
                     column_title = main.title,
                     column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                     heatmap_legend_side = "right",
                     padding = grid::unit(c(2, 2, 2, 2), "mm"))

    ## Print interpretation (if verbose)
    if (verbose) {
        cat("Interpretation:\n")
        if (!is.null(cluster.stats)) {
            if (cluster.stats$rank1.variance > 0.90 && cluster.stats$effective.rank == 1) {
                cat("  -> Simple structure: rank-1 model appropriate\n")
                cat("  -> Phylotypes show coherent positive/negative pattern\n")
            } else if (cluster.stats$effective.rank <= 3) {
                cat("  -> Moderate complexity: rank-2 or rank-3 model needed\n")
                cat("  -> Multiple independent co-monotonicity patterns present\n")
            } else {
                cat("  -> High complexity: rank-", cluster.stats$effective.rank, "model needed\n")
                cat("  -> Multiple independent patterns - NOT suitable for rank-1!\n")
                cat("  -> This is where BCs4vd would cherry-pick phylotypes\n")
            }
        }
        cat("\n")
    }


    ## Return results
    return(invisible(list(
        heatmap = ht,
        drawn = ht.drawn,
        cluster.id = cluster.id,
        cluster.data = cluster.data,
        cluster.members = cluster.members,
        cluster.info = cluster.info,
        cluster.stats = cluster.stats,
        phylo.means = phylo.means,
        phylo.sds = phylo.sds,
        strong.positive = names(phylo.means[phylo.means > phylo.mean.thld2]),
        strong.negative = names(phylo.means[phylo.means < -phylo.mean.thld2]),
        n.phylotypes.original = n.phylotypes.original,
        n.phylotypes.displayed = ncol(cluster.data),
        phylo.filter.thld = phylo.filter.thld,
        phylo.mean.thld1 = phylo.mean.thld1,
        phylo.mean.thld2 = phylo.mean.thld2
    )))
}


#' Compare Multiple MGCP Clusters Side-by-Side
#'
#' Creates a multi-panel visualization comparing several MGCP clusters.
#'
#' @param data.matrix Numeric matrix (n.vertices Ã— n.phylotypes)
#' @param mgcp.result Result from \code{extract.stable.clusters()}
#' @param cluster.ids Integer vector specifying which clusters to compare
#' @param cluster.rows Logical for hierarchical clustering within clusters
#' @param cluster.cols Logical for hierarchical clustering of phylotypes
#' @param show.phylotype.names Logical for phylotype name display
#' @param ncol.layout Number of columns in layout (default: 2)
#' @param verbose Logical indicating whether to print progress messages. Default is \code{TRUE}.
#'
#' @return List of heatmap objects
#'
#' @examples
#' \dontrun{
#' # Compare top 4 most stable clusters
#' compare.mgcp.clusters(comono.smooth, mgcp.result, cluster.ids = 1:4)
#'
#' # Compare clusters with different complexity
#' compare.mgcp.clusters(comono.smooth, mgcp.result, cluster.ids = c(1, 4, 12))
#' }
#'
#' @export
compare.mgcp.clusters <- function(data.matrix,
                                  mgcp.result,
                                  cluster.ids,
                                  cluster.rows = TRUE,
                                  cluster.cols = TRUE,
                                  show.phylotype.names = FALSE,
                                  ncol.layout = 2,
                                  verbose = TRUE) {

  if (length(cluster.ids) < 2) {
    stop("Need at least 2 clusters to compare")
  }

  if (length(cluster.ids) > 6) {
    warning("Comparing > 6 clusters may create crowded visualization")
  }

  if (verbose) {
    cat("\n")
    cat("Comparing", length(cluster.ids), "MGCP clusters\n")
    cat(strrep("=", 50), "\n\n")
  }

  # Create list of heatmaps
  ht.list <- list()

  for (i in seq_along(cluster.ids)) {
    cl.id <- cluster.ids[i]

    if (verbose) {
      cat("Creating visualization for Cluster", cl.id, "...\n")
    }

    # Capture the heatmap without drawing
    ht.result <- visualize.mgcp.cluster(
      data.matrix,
      mgcp.result,
      cluster.id = cl.id,
      cluster.rows = cluster.rows,
      cluster.cols = cluster.cols,
      show.phylotype.names = show.phylotype.names,
      show.statistics = FALSE,  # Too crowded in comparison
      verbose = FALSE  # Suppress verbose output for individual clusters
    )

    ht.list[[i]] <- ht.result$heatmap
  }

  if (verbose) {
    cat("\n")
  }

  # Determine layout
  n.clusters <- length(cluster.ids)
  nrow.layout <- ceiling(n.clusters / ncol.layout)

  # Create combined plot
  if (verbose) {
    cat("Drawing combined visualization...\n")
  }

  # Use grid to create layout
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow.layout, ncol.layout)))

  for (i in seq_along(ht.list)) {
    row.idx <- ceiling(i / ncol.layout)
    col.idx <- ((i - 1) %% ncol.layout) + 1

    grid::pushViewport(grid::viewport(layout.pos.row = row.idx,
                         layout.pos.col = col.idx))
    ComplexHeatmap::draw(ht.list[[i]], newpage = FALSE)
    grid::upViewport()
  }

  grid::upViewport()

  if (verbose) {
    cat("Comparison complete!\n\n")

    # Print summary comparison
    cat("Cluster Comparison Summary:\n")
    cat(strrep("-", 50), "\n")

    cluster.analysis <- analyze.extracted.clusters(data.matrix, mgcp.result)

    for (cl.id in cluster.ids) {
      cl.info <- mgcp.result$summary.stats[mgcp.result$summary.stats$cluster.id == cl.id, ]
      cl.stats <- cluster.analysis[cluster.analysis$cluster.id == cl.id, ]

      cat(sprintf("Cluster %d: n=%d, gap=%.1f, rank-1=%.1f%%, eff.rank=%d\n",
                  cl.id,
                  cl.info$cluster.size,
                  cl.info$gap.size,
                  cl.stats$rank1.variance * 100,
                  cl.stats$effective.rank))
    }
    cat("\n")
  }

  return(invisible(ht.list))
}


#' Export MGCP Cluster Data
#'
#' Exports a single cluster's data and statistics to CSV files for further analysis.
#'
#' @param data.matrix Original data matrix
#' @param mgcp.result Result from \code{extract.stable.clusters()}
#' @param cluster.id Cluster to export
#' @param output.prefix Prefix for output filenames
#' @param verbose Logical indicating whether to print export messages. Default is \code{TRUE}.
#' @param phylo.mean.thld2 Numeric threshold for strong phylotype signal in exported
#'   statistics file. Default is 0.3.
#'
#' @return Invisibly returns list of exported data
#'
#' @examples
#' \dontrun{
#' # Export cluster 1 data
#' export.mgcp.cluster(comono.smooth, mgcp.result, cluster.id = 1,
#'                    output.prefix = "cluster1")
#' # Creates: cluster1_data.csv, cluster1_stats.csv, cluster1_phylotypes.csv
#' }
#'
#' @export
export.mgcp.cluster <- function(data.matrix,
                                mgcp.result,
                                cluster.id,
                                output.prefix = NULL,
                                verbose = TRUE,
                                phylo.mean.thld2 = 0.3) {

  if (is.null(output.prefix)) {
    output.prefix <- sprintf("mgcp_cluster_%d", cluster.id)
  }

  # Extract cluster data
  cluster.members <- which(mgcp.result$clusters == cluster.id)
  cluster.data <- data.matrix[cluster.members, , drop = FALSE]

  # Get statistics
  cluster.info <- mgcp.result$summary.stats[mgcp.result$summary.stats$cluster.id == cluster.id, ]
  cluster.analysis <- analyze.extracted.clusters(data.matrix, mgcp.result)
  cluster.stats <- cluster.analysis[cluster.analysis$cluster.id == cluster.id, ]

  # Phylotype statistics
  phylo.means <- colMeans(cluster.data)
  phylo.sds <- apply(cluster.data, 2, sd)
  phylo.stats <- data.frame(
    phylotype = colnames(cluster.data),
    mean.value = phylo.means,
    sd.value = phylo.sds,
    min.value = apply(cluster.data, 2, min),
    max.value = apply(cluster.data, 2, max),
    strong.signal = abs(phylo.means) > phylo.mean.thld2,
    stringsAsFactors = FALSE
  )

  # Export files
  write.csv(cluster.data,
           paste0(output.prefix, "_data.csv"),
           row.names = TRUE)

  write.csv(cbind(cluster.info, cluster.stats),
           paste0(output.prefix, "_stats.csv"),
           row.names = FALSE)

  write.csv(phylo.stats,
           paste0(output.prefix, "_phylotypes.csv"),
           row.names = FALSE)

  if (verbose) {
    cat("Exported cluster", cluster.id, "data:\n")
    cat("  ", paste0(output.prefix, "_data.csv"), "\n")
    cat("  ", paste0(output.prefix, "_stats.csv"), "\n")
    cat("  ", paste0(output.prefix, "_phylotypes.csv"), "\n\n")
  }

  return(invisible(list(
    data = cluster.data,
    info = cluster.info,
    stats = cluster.stats,
    phylotypes = phylo.stats
  )))
}

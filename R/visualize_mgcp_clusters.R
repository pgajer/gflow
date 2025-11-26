#' Visualize Maximal Gap Cluster Peeling Results
#'
#' Creates a comprehensive heatmap visualization of MGCP clustering results.
#' Unlike bi-clustering methods, MGCP produces genuine disjoint vertex clusters
#' without feature filtering, so all phylotypes are displayed.
#'
#' @param data.matrix Numeric matrix to be visualized (n.vertices × n.phylotypes).
#' @param mgcp.result Result object from \code{extract.stable.clusters()} containing
#'   MGCP clustering results.
#' @param cluster.rows Logical indicating whether to perform hierarchical clustering
#'   on rows within clusters. Default is \code{FALSE} to preserve extraction order.
#'   When \code{TRUE}, rows are clustered **within each cluster block** (not globally),
#'   which can reveal internal structure while maintaining cluster separation.
#' @param cluster.cols Logical indicating whether to perform hierarchical clustering
#'   on columns. Default is \code{TRUE}.
#' @param show.unassigned Logical indicating whether to display unassigned vertices
#'   (cluster ID = 0). Default is \code{TRUE}.
#' @param show.cluster.labels Logical indicating whether to add cluster numbers
#'   directly next to the vertex cluster color bars. Default is \code{TRUE}.
#' @param show.extraction.order Logical indicating whether to sort clusters by
#'   extraction order (most stable first). If \code{FALSE}, sorts by cluster size.
#'   Default is \code{TRUE}.
#' @param color.palette Function or character vector for cluster colors. Default
#'   uses rainbow colors.
#' @param unassigned.color Color for unassigned vertices. Default is "gray80".
#' @param main.title Character string for heatmap title. If NULL, generates
#'   automatic title with cluster count.
#' @param show.phylotype.names Logical indicating whether to display phylotype
#'   (feature) names at the bottom of the heatmap. Default is \code{TRUE}.
#'   Set to \code{FALSE} for cleaner display when many phylotypes present.
#' @param row.title Character string for row axis title. Default is "Vertices".
#' @param column.title Character string for column axis title. Default is "Phylotypes".
#' @param heatmap.legend.title Heatmap legend title.
#' @param verbose Logical indicating whether to print summary statistics and
#'   cluster information to console. Default is \code{TRUE}.
#'
#' @return Returns invisibly a list containing heatmap object, ordered data, and
#'   cluster information.
#'
#' @details
#' This function visualizes MGCP clustering results with several key features:
#' \itemize{
#'   \item All phylotypes (features) are retained - no filtering
#'   \item Clusters are genuinely disjoint - no overlap handling needed
#'   \item Extraction order reflects cluster stability (gap size)
#'   \item Unassigned vertices can be shown or hidden
#'   \item Cluster labels displayed inline for easy identification
#' }
#'
#' The default ordering by extraction order is meaningful: Cluster 1 is the most
#' stable (largest gap), Cluster 2 is second most stable, etc. This provides
#' immediate insight into cluster quality.
#'
#' @seealso \code{\link{extract.stable.clusters}},
#'   \code{\link{analyze.extracted.clusters}}
#'
#' @examples
#' \dontrun{
#' # Extract clusters using MGCP
#' mgcp.result <- extract.stable.clusters(
#'   comono.smooth,
#'   min.cluster.size = 20
#' )
#'
#' # Visualize with default settings (phylotype names shown)
#' ht <- visualize.mgcp.clusters(comono.smooth, mgcp.result)
#'
#' # Hide phylotype names for cleaner display
#' ht <- visualize.mgcp.clusters(
#'   comono.smooth,
#'   mgcp.result,
#'   show.phylotype.names = FALSE
#' )
#'
#' # Hide unassigned vertices, sort by size, hide phylotype names, silent operation
#' ht <- visualize.mgcp.clusters(
#'   comono.smooth,
#'   mgcp.result,
#'   show.unassigned = FALSE,
#'   show.extraction.order = FALSE,
#'   show.phylotype.names = FALSE,
#'   verbose = FALSE
#' )
#'
#' # Custom axis titles
#' ht <- visualize.mgcp.clusters(
#'   comono.smooth,
#'   mgcp.result,
#'   row.title = "Samples",
#'   column.title = "Features (all retained)"
#' )
#' }
#'
#' @export
visualize.mgcp.clusters <- function(data.matrix,
                                    mgcp.result,
                                    cluster.rows = FALSE,
                                    cluster.cols = TRUE,
                                    show.unassigned = TRUE,
                                    show.cluster.labels = TRUE,
                                    show.extraction.order = TRUE,
                                    color.palette = NULL,
                                    unassigned.color = "gray80",
                                    main.title = NULL,
                                    show.phylotype.names = TRUE,
                                    row.title = "Vertices",
                                    column.title = "Phylotypes",
                                    heatmap.legend.title = "Co-monotonicity coefficients",
                                    verbose = TRUE) {

  # Extract cluster assignments
  cluster.assignments <- mgcp.result$clusters
  n.clusters <- mgcp.result$n.clusters
  n.total <- length(cluster.assignments)

  if (n.clusters == 0) {
    stop("No clusters found in mgcp.result")
  }

  # Identify assigned and unassigned vertices
  assigned.vertices <- cluster.assignments > 0
  n.assigned <- sum(assigned.vertices)
  n.unassigned <- sum(!assigned.vertices)

  if (verbose) {
    cat("MGCP Visualization Summary:\n")
    cat("  Total vertices:", n.total, "\n")
    cat("  Clusters extracted:", n.clusters, "\n")
    cat("  Vertices assigned:", n.assigned, "\n")
    cat("  Vertices unassigned:", n.unassigned, "\n")
    cat("  Total phylotypes:", ncol(data.matrix), "(all retained)\n\n")
  }

  # Determine which vertices to display
  if (show.unassigned) {
    display.vertices <- 1:n.total
    display.assignments <- cluster.assignments
  } else {
    display.vertices <- which(assigned.vertices)
    display.assignments <- cluster.assignments[assigned.vertices]

    if (verbose) {
      cat("Excluding unassigned vertices from display\n")
      cat("  Displaying:", length(display.vertices), "vertices\n\n")
    }
  }

  if (length(display.vertices) == 0) {
    stop("No vertices to display")
  }

  # Prepare cluster ordering
  cluster.ids <- sort(unique(display.assignments[display.assignments > 0]))
  n.display.clusters <- length(cluster.ids)

  # Get cluster information for ordering/labeling
  cluster.info <- mgcp.result$summary.stats
  cluster.info <- cluster.info[cluster.info$cluster.id %in% cluster.ids, ]

  if (show.extraction.order) {
    # Order by extraction order (already sorted in summary.stats)
    cluster.order <- cluster.info$cluster.id
    ordering.type <- "extraction order (stability)"
  } else {
    # Order by size
    cluster.order <- cluster.info$cluster.id[order(-cluster.info$cluster.size)]
    ordering.type <- "cluster size"
  }

  if (verbose) {
    cat("Cluster ordering:", ordering.type, "\n")
    for (i in seq_along(cluster.order)) {
      cl.id <- cluster.order[i]
      cl.info <- cluster.info[cluster.info$cluster.id == cl.id, ]
      cat("  Position", i, ": Cluster", cl.id,
          "- size =", cl.info$cluster.size,
          ", gap =", round(cl.info$gap.size, 2), "\n")
    }
    cat("\n")
  }

  # Create ordered data matrix
  # Sort vertices by cluster membership according to chosen ordering
  vertex.order.factor <- factor(display.assignments,
                                levels = c(cluster.order, if(show.unassigned) 0 else NULL))
  vertex.order <- order(vertex.order.factor)

  ordered.matrix <- data.matrix[display.vertices[vertex.order], , drop = FALSE]
  ordered.assignments <- display.assignments[vertex.order]

  # Set up colors
  if (is.null(color.palette)) {
    cluster.colors <- rainbow(n.display.clusters)
  } else if (is.function(color.palette)) {
    cluster.colors <- color.palette(n.display.clusters)
  } else {
    cluster.colors <- rep_len(color.palette, n.display.clusters)
  }

  # Create color mapping
  color.mapping <- setNames(cluster.colors, as.character(cluster.order))
  if (show.unassigned) {
    color.mapping <- c(color.mapping, "0" = unassigned.color)
  }

  # Prepare cluster labels
  cluster.labels <- character(length(unique(ordered.assignments)))
  names(cluster.labels) <- as.character(sort(unique(ordered.assignments)))

  for (cl in cluster.order) {
    cl.info <- cluster.info[cluster.info$cluster.id == cl, ]
    cluster.labels[as.character(cl)] <- sprintf(
      "C%d (n=%d, gap=%.1f)",
      cl,
      cl.info$cluster.size,
      cl.info$gap.size
    )
  }

  if (show.unassigned) {
    cluster.labels["0"] <- sprintf("Unassigned (n=%d)", n.unassigned)
  }

  # Create row annotation
  row.ha <- ComplexHeatmap::rowAnnotation(
    Cluster = factor(ordered.assignments,
                    levels = if(show.unassigned) c(cluster.order, 0) else cluster.order),
    col = list(Cluster = color.mapping),
    show_legend = FALSE,  # We'll use inline labels instead
    annotation_name_side = "top",
    simple_anno_size = grid::unit(6, "mm")
  )

  # Color scheme for co-monotonicity
  col.fun <- circlize::colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    c("blue", "lightblue", "white", "pink", "red")
  )

  # Generate title
  if (is.null(main.title)) {
    main.title <- sprintf(
      "MGCP: %d Clusters (ordered by %s)",
      n.clusters,
      if(show.extraction.order) "stability" else "size"
    )
  }

  # Create heatmap
  # IMPORTANT: cluster_row_slices = FALSE ensures that when cluster.rows = TRUE,
  # hierarchical clustering happens WITHIN each cluster block, not globally.
  # Without this, ComplexHeatmap would reorder the cluster blocks themselves,
  # disrupting our intended ordering (by stability or size).
  ht <- ComplexHeatmap::Heatmap(
    ordered.matrix,
    name = "Co-monotonicity",
    col = col.fun,
    cluster_rows = cluster.rows,
    cluster_columns = cluster.cols,
    cluster_row_slices = FALSE,  # Keep cluster blocks in specified order
    cluster_column_slices = FALSE,  # Don't reorder column splits if any
    show_row_names = FALSE,
    show_column_names = show.phylotype.names,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_rot = 45,
    column_names_side = "bottom",
    row_title = row.title,
    row_title_side = "left",
    column_title = column.title,
    column_title_side = "top",
    right_annotation = row.ha,
    row_split = factor(ordered.assignments,
                      levels = if(show.unassigned) c(cluster.order, 0) else cluster.order),
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    row_gap = grid::unit(2, "mm"),
    border = TRUE,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = heatmap.legend.title,
      direction = "horizontal",
      legend_width = grid::unit(6, "cm")
    ),
    use_raster = nrow(ordered.matrix) > 1000 || ncol(ordered.matrix) > 1000,
    raster_quality = 2
  )

  # Add cluster labels inline
  if (show.cluster.labels) {
    label.levels <- if(show.unassigned) c(cluster.order, 0) else cluster.order
    ht <- ht + ComplexHeatmap::rowAnnotation(
      Labels = ComplexHeatmap::anno_block(
        gp = grid::gpar(fill = NA, col = NA),
        labels = cluster.labels[as.character(label.levels)],
        labels_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        labels_rot = 0,
        which = "row"
      ),
      width = grid::unit(35, "mm")
    )
  }

  # Draw heatmap
  ht.drawn <- ComplexHeatmap::draw(ht, 
                   column_title = main.title,
                   heatmap_legend_side = "bottom",
                   padding = grid::unit(c(2, 2, 2, 10), "mm"))

  # Print cluster statistics
  if (verbose) {
    cat("\n")
    cat("Cluster Statistics:\n")
    cat("==================\n")
    print(cluster.info[, c("cluster.id", "cluster.size", "cluster.height", "gap.size")])
  }

  return(invisible(list(
    heatmap = ht,
    drawn = ht.drawn,
    ordered.matrix = ordered.matrix,
    ordered.assignments = ordered.assignments,
    vertex.order = display.vertices[vertex.order],
    cluster.info = cluster.info,
    color.mapping = color.mapping,
    n.clusters = n.clusters,
    n.unassigned = n.unassigned
  )))
}

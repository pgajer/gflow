#' Create Heatmap with Labeled Super-Cell Partition Color Bar
#'
#' Generates a heatmap using ComplexHeatmap with cell numbers displayed in the
#' internal margin between the color bar and heatmap body. All super-cells are
#' guaranteed to have labels.
#'
#' @param cm Numeric matrix (co-monotonicity or similar). Rows represent samples.
#' @param partition Numeric vector of super-cell membership for rows, same length
#'   as nrow(cm).
#' @param color.palette A circlize colorRamp2 object or NULL for default blue-yellow-red.
#' @param cluster.cols Logical indicating whether to cluster columns. Default is TRUE.
#' @param main Character string for plot title.
#' @param show.row.names Logical indicating whether to show row names. Default is FALSE.
#' @param show.col.names Logical indicating whether to show column names. Default is FALSE.
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap().
#'
#' @return A ComplexHeatmap object.
#'
#' @details
#' The function organizes rows according to their super-cell partition membership,
#' with hierarchical clustering (Ward's method) applied within each super-cell.
#' Cell numbers appear in the margin between the annotation bar and heatmap,
#' ensuring all cells are clearly labeled.
#'
#' @examples
#' \dontrun{
#' library(circlize)
#'
#' ## Create color palette
#' blue.yellow.red.palette <- colorRamp2(c(-1, 0, 1), c("blue", "yellow", "red"))
#'
#' ## Generate heatmap
#' ht <- partition.heatmap(
#'   cm = sptb.spp.cm,
#'   partition = sptb.phylo.super.cells,
#'   color.palette = blue.yellow.red.palette,
#'   main = "sPTB Species Co-Monotonicity Matrix"
#' )
#'
#' ## Display
#' draw(ht)
#'
#' ## Save to PDF
#' pdf("heatmap.pdf", width = 8, height = 10)
#' draw(ht)
#' dev.off()
#' }
#'
#' @export
partition.heatmap <- function(cm,
                              partition,
                              color.palette = NULL,
                              cluster.cols = TRUE,
                              main = "Co-Monotonicity Matrix",
                              show.row.names = FALSE,
                              show.col.names = FALSE,
                              ...) {

  ## Load required packages
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Install from Bioconductor:\n",
         "  BiocManager::install('ComplexHeatmap')")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Install with:\n",
         "  install.packages('circlize')")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required")
  }

  ## Input validation
  if (!is.matrix(cm) && !is.data.frame(cm)) {
    stop("cm must be a matrix or data frame")
  }
  if (!is.numeric(partition)) {
    stop("partition must be a numeric vector")
  }
  if (length(partition) != nrow(cm)) {
    stop("partition length must equal number of rows in cm")
  }
  if (any(is.na(partition))) {
    stop("partition cannot contain NA values")
  }

  ## Convert to matrix if needed
  if (is.data.frame(cm)) {
    cm <- as.matrix(cm)
  }

  ## Default color palette if not provided
  if (is.null(color.palette)) {
    color.palette <- circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "yellow", "red")
    )
  }

  ## Get unique super-cell IDs
  cell.ids <- sort(unique(partition))
  n.cells <- length(cell.ids)

  ## Compute row ordering: cluster within each super-cell
  row.order <- integer(0)

  for (cell.id in cell.ids) {
    cell.rows <- which(partition == cell.id)

    if (length(cell.rows) == 1) {
      row.order <- c(row.order, cell.rows)
    } else {
      cell.matrix <- cm[cell.rows, , drop = FALSE]
      row.dist <- dist(cell.matrix)
      row.hclust <- hclust(row.dist, method = "ward.D2")
      row.order <- c(row.order, cell.rows[row.hclust$order])
    }
  }

  ## Reorder matrix
  cm.ordered <- cm[row.order, , drop = FALSE]
  partition.ordered <- partition[row.order]

  ## Create color palette for super-cells
  if (n.cells <= 12) {
    cell.colors <- RColorBrewer::brewer.pal(max(3, n.cells), "Set3")[1:n.cells]
  } else {
    cell.colors <- rainbow(n.cells)
  }
  names(cell.colors) <- as.character(cell.ids)

  ## Create row annotation
  row.annotation <- ComplexHeatmap::rowAnnotation(
    Cell = partition.ordered,
    col = list(Cell = cell.colors),
    annotation_label = "Super-Cell",
    annotation_name_side = "top",
    show_legend = FALSE,
    simple_anno_size = grid::unit(0.5, "cm"),
    annotation_name_gp = grid::gpar(fontsize = 10)
  )

  ## Create column clustering if requested
  if (cluster.cols) {
    col.cluster <- TRUE
  } else {
    col.cluster <- FALSE
  }

  ## Create row titles (one per unique cell ID, in order)
  ## This must match the number of splits created by row_split
  row.titles <- as.character(cell.ids)

  ## Create the heatmap
  ht <- ComplexHeatmap::Heatmap(
    cm.ordered,
    name = "Value",
    col = color.palette,
    cluster_rows = FALSE,
    cluster_columns = col.cluster,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    show_row_names = show.row.names,
    show_column_names = show.col.names,
    column_title = main,
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    left_annotation = row.annotation,
    row_split = factor(partition.ordered, levels = cell.ids),  ## Factor with ordered levels
    row_title = row.titles,  ## One title per split level
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    row_gap = grid::unit(2, "mm"),
    border = TRUE,
    use_raster = TRUE,  ## Explicitly set for large matrices
    raster_quality = 2,  ## Higher quality rasterization
    heatmap_legend_param = list(
      title = "Value",
      at = c(-1, 0, 1),
      labels = c("-1", "0", "1"),
      legend_height = grid::unit(4, "cm"),
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9)
    ),
    ...
  )

  return(ht)
}

#' Create Individual Heatmaps for Each Super-Cell
#'
#' Generates separate heatmaps for each super-cell in a partition, showing
#' hierarchical clustering of both samples (rows) and phylotypes (columns)
#' within each cell.
#'
#' @param cm Numeric matrix (co-monotonicity or similar). Rows represent samples,
#'   columns represent phylotypes.
#' @param partition Numeric vector of super-cell membership for rows, same length
#'   as nrow(cm).
#' @param color.palette A circlize colorRamp2 object or NULL for default blue-yellow-red.
#' @param cell.ids Optional vector of specific cell IDs to plot. If NULL, plots all cells.
#' @param show.row.names Logical indicating whether to show row names. Default is FALSE.
#' @param show.col.names Logical indicating whether to show column names. Default is FALSE.
#' @param main.prefix Character string prefix for plot titles. Cell ID will be appended.
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap().
#'
#' @return A named list of ComplexHeatmap objects, one per super-cell.
#'
#' @details
#' The function creates individual heatmaps for each super-cell, allowing detailed
#' examination of the co-monotonicity structure within each cell. Both rows (samples)
#' and columns (phylotypes) are hierarchically clustered using Ward's method.
#'
#' @examples
#' \dontrun{
#' library(ComplexHeatmap)
#' library(circlize)
#'
#' ## Create color palette
#' blue.yellow.red.palette <- colorRamp2(c(-1, 0, 1), c("blue", "yellow", "red"))
#'
#' ## Generate heatmaps for all cells
#' cell.heatmaps <- partition.cell.heatmaps(
#'   cm = sptb.spp.cm,
#'   partition = sptb.phylo.super.cells,
#'   color.palette = blue.yellow.red.palette
#' )
#'
#' ## Display a specific cell
#' draw(cell.heatmaps[["3"]])
#'
#' ## Save all to PDF
#' pdf("cell_heatmaps.pdf", width = 10, height = 8)
#' for (cell.id in names(cell.heatmaps)) {
#'   draw(cell.heatmaps[[cell.id]])
#' }
#' dev.off()
#'
#' ## Or save each cell to separate PDF
#' for (cell.id in names(cell.heatmaps)) {
#'   pdf(paste0("cell_", cell.id, "_heatmap.pdf"), width = 10, height = 8)
#'   draw(cell.heatmaps[[cell.id]])
#'   dev.off()
#' }
#' }
#'
#' @export
partition.cell.heatmaps <- function(cm,
                                    partition,
                                    color.palette = NULL,
                                    cell.ids = NULL,
                                    show.row.names = FALSE,
                                    show.col.names = FALSE,
                                    main.prefix = "Super-Cell",
                                    ...) {

  ## Load required packages
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Install from Bioconductor:\n",
         "  BiocManager::install('ComplexHeatmap')")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Install with:\n",
         "  install.packages('circlize')")
  }

  ## Input validation
  if (!is.matrix(cm) && !is.data.frame(cm)) {
    stop("cm must be a matrix or data frame")
  }
  if (!is.numeric(partition)) {
    stop("partition must be a numeric vector")
  }
  if (length(partition) != nrow(cm)) {
    stop("partition length must equal number of rows in cm")
  }
  if (any(is.na(partition))) {
    stop("partition cannot contain NA values")
  }

  ## Convert to matrix if needed
  if (is.data.frame(cm)) {
    cm <- as.matrix(cm)
  }

  ## Default color palette if not provided
  if (is.null(color.palette)) {
    color.palette <- circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "yellow", "red")
    )
  }

  ## Get cell IDs to plot
  if (is.null(cell.ids)) {
    cell.ids <- sort(unique(partition))
  } else {
    ## Validate requested cell IDs
    available.cells <- unique(partition)
    invalid.cells <- setdiff(cell.ids, available.cells)
    if (length(invalid.cells) > 0) {
      stop("Requested cell IDs not found in partition: ",
           paste(invalid.cells, collapse = ", "))
    }
  }

  ## Create list to store heatmaps
  heatmap.list <- list()

  ## Generate heatmap for each cell
  for (cell.id in cell.ids) {
    ## Get rows for this cell
    cell.rows <- which(partition == cell.id)
    n.rows <- length(cell.rows)

    ## Subset matrix
    cell.matrix <- cm[cell.rows, , drop = FALSE]

    ## Determine if we should use raster (for large matrices)
    use.raster <- n.rows > 500

    ## Create title
    cell.title <- paste0(main.prefix, " ", cell.id, " (n = ", n.rows, ")")

    ## Create heatmap
    ht <- ComplexHeatmap::Heatmap(
      cell.matrix,
      name = "Value",
      col = color.palette,
      cluster_rows = TRUE,  ## Cluster samples within cell
      cluster_columns = TRUE,  ## Cluster phylotypes
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "ward.D2",
      clustering_method_columns = "ward.D2",
      show_row_names = show.row.names,
      show_column_names = show.col.names,
      column_title = cell.title,
      column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      row_title = NULL,
      column_names_gp = grid::gpar(fontsize = 8),
      row_names_gp = grid::gpar(fontsize = 8),
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      row_dend_width = grid::unit(2, "cm"),
      column_dend_height = grid::unit(2, "cm"),
      border = TRUE,
      use_raster = use.raster,
      raster_quality = 2,
      heatmap_legend_param = list(
        title = "Value",
        at = c(-1, 0, 1),
        labels = c("-1", "0", "1"),
        legend_height = grid::unit(4, "cm"),
        title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9)
      ),
      ...
    )

    ## Store in list
    heatmap.list[[as.character(cell.id)]] <- ht
  }

  return(heatmap.list)
}

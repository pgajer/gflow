#' Visualize Bi-clustering Results with Optional Filtering and Disjoint Assignment
#'
#' Creates a comprehensive heatmap visualization of bi-clustering results with options
#' to create disjoint vertex assignments, filter phylotypes, and control annotation display.
#' This produces a cleaner high-level overview of co-monotonicity patterns.
#'
#' @param data.matrix Numeric matrix to be visualized (n.vertices x n.phylotypes).
#' @param biclust.result An object of class \code{Biclust} containing bi-clustering results.
#' @param max.biclusters Integer specifying maximum number of bi-clusters to display.
#'   Default is 10.
#' @param cluster.rows Logical indicating whether to perform hierarchical clustering
#'   on rows within bi-clusters. Default is \code{TRUE}.
#' @param cluster.cols Logical indicating whether to perform hierarchical clustering
#'   on columns within bi-clusters. Default is \code{TRUE}.
#' @param make.disjoint Logical indicating whether to make vertex (row) assignments
#'   disjoint. Default is \code{FALSE}.
#' @param filter.phylotypes Logical indicating whether to remove phylotypes that do
#'   not appear in any bi-cluster. Default is \code{FALSE}.
#' @param assignment.priority Character string specifying assignment priority when
#'   \code{make.disjoint = TRUE}. Options: "first", "largest", "smallest".
#'   Default is "first".
#' @param show.phylotype.biclusters Logical indicating whether to show phylotype
#'   bi-cluster membership as colored bars at top of heatmap. If \code{FALSE}, only
#'   the phylotype dendrogram is shown. Default is \code{TRUE}.
#' @param show.bicluster.labels Logical indicating whether to add bi-cluster numbers
#'   directly next to the vertex bi-cluster color bars instead of using a separate
#'   legend. Default is \code{FALSE}.
#'
#' @return Returns invisibly a list containing heatmap object, filtered data, and
#'   assignment information.
#'
#' @details
#' This enhanced visualization function provides multiple options for creating cleaner
#' overview heatmaps. Setting \code{show.phylotype.biclusters = FALSE} removes the
#' top annotation bars showing phylotype membership, creating a simpler view focused
#' on vertex patterns. Setting \code{show.bicluster.labels = TRUE} adds bi-cluster
#' numbers directly next to the colored bars on the right, eliminating the need for
#' a separate legend and making it easier to identify specific bi-clusters.
#'
#' @seealso \code{\link{visualize.bcs4vd.bicluster}},
#'   \code{\link{summarize.bcs4vd.biclusters}}
#'
#' @examples
#' \dontrun{
#' # Clean visualization with phylotype bars hidden and inline labels
#' ht <- visualize.bcs4vd.biclusters(
#'   comono.smooth,
#'   bc.result,
#'   max.biclusters = 10,
#'   make.disjoint = TRUE,
#'   filter.phylotypes = TRUE,
#'   show.phylotype.biclusters = FALSE,
#'   show.bicluster.labels = TRUE
#' )
#' }
#'
#' @export
visualize.bcs4vd.biclusters <- function(data.matrix,
                                         biclust.result,
                                         max.biclusters = 10,
                                         cluster.rows = TRUE,
                                         cluster.cols = TRUE,
                                         make.disjoint = FALSE,
                                         filter.phylotypes = FALSE,
                                         assignment.priority = c("first", "largest", "smallest"),
                                         show.phylotype.biclusters = TRUE,
                                         show.bicluster.labels = FALSE) {

  assignment.priority <- match.arg(assignment.priority)
  n.bc <- min(biclust.result@Number, max.biclusters)

  if (n.bc == 0) {
    stop("No bi-clusters found in biclust.result")
  }

  # Filter phylotypes if requested
  phylotype.in.any.bc <- colSums(biclust.result@NumberxCol) > 0

  if (filter.phylotypes) {
    phylotype.indices <- which(phylotype.in.any.bc)
    n.phylotypes.removed <- sum(!phylotype.in.any.bc)

    cat("Filtering phylotypes:\n")
    cat("  Total phylotypes:", ncol(data.matrix), "\n")
    cat("  Phylotypes in bi-clusters:", length(phylotype.indices), "\n")
    cat("  Phylotypes removed:", n.phylotypes.removed, "\n\n")

    if (length(phylotype.indices) == 0) {
      stop("No phylotypes found in any bi-cluster")
    }
  } else {
    phylotype.indices <- 1:ncol(data.matrix)
    n.phylotypes.removed <- 0
  }

  # Handle vertex assignment (disjoint or overlapping)
  n.vertices <- nrow(data.matrix)
  vertex.assignment <- rep(0, n.vertices)

  if (make.disjoint) {
    cat("Creating disjoint vertex assignments:\n")

    bc.sizes <- colSums(biclust.result@RowxNumber[, 1:n.bc])

    if (assignment.priority == "first") {
      for (i in 1:n.bc) {
        vertices.in.bc <- which(biclust.result@RowxNumber[, i])
        unassigned <- vertices.in.bc[vertex.assignment[vertices.in.bc] == 0]
        vertex.assignment[unassigned] <- i
      }

    } else if (assignment.priority == "largest") {
      bc.order <- order(bc.sizes, decreasing = TRUE)
      for (i in bc.order) {
        vertices.in.bc <- which(biclust.result@RowxNumber[, i])
        unassigned <- vertices.in.bc[vertex.assignment[vertices.in.bc] == 0]
        vertex.assignment[unassigned] <- i
      }

    } else if (assignment.priority == "smallest") {
      bc.order <- order(bc.sizes, decreasing = FALSE)
      for (i in bc.order) {
        vertices.in.bc <- which(biclust.result@RowxNumber[, i])
        unassigned <- vertices.in.bc[vertex.assignment[vertices.in.bc] == 0]
        vertex.assignment[unassigned] <- i
      }
    }

    n.assigned <- sum(vertex.assignment > 0)
    n.vertices.removed <- n.vertices - n.assigned

    cat("  Assignment priority:", assignment.priority, "\n")
    cat("  Vertices assigned:", n.assigned, "\n")
    cat("  Vertices not in any bi-cluster:", n.vertices.removed, "\n")

    bc.counts <- table(vertex.assignment[vertex.assignment > 0])
    cat("  Vertices per bi-cluster (after disjoint assignment):\n")
    for (bc.id in names(bc.counts)) {
      orig.count <- sum(biclust.result@RowxNumber[, as.integer(bc.id)])
      cat("    BC", bc.id, ":", bc.counts[bc.id],
          "(was", orig.count, ", removed", orig.count - bc.counts[bc.id], ")\n")
    }
    cat("\n")

    vertex.indices <- which(vertex.assignment > 0)

  } else {
    vertex.in.any.bc <- rowSums(biclust.result@RowxNumber[, 1:n.bc]) > 0
    vertex.indices <- which(vertex.in.any.bc)
    n.vertices.removed <- sum(!vertex.in.any.bc)

    for (i in 1:n.bc) {
      vertices.in.bc <- which(biclust.result@RowxNumber[, i])
      unassigned <- vertices.in.bc[vertex.assignment[vertices.in.bc] == 0]
      vertex.assignment[unassigned] <- i
    }

    cat("Keeping overlapping vertex assignments:\n")
    cat("  Vertices in bi-clusters:", length(vertex.indices), "\n")
    cat("  Vertices removed:", n.vertices.removed, "\n\n")
  }

  # Filter matrix
  filtered.matrix <- data.matrix[vertex.indices, phylotype.indices, drop = FALSE]
  vertex.assignment.filtered <- vertex.assignment[vertex.indices]

  # Order vertices by bi-cluster membership
  vertex.order <- order(vertex.assignment.filtered)
  filtered.matrix <- filtered.matrix[vertex.order, , drop = FALSE]
  vertex.assignment.ordered <- vertex.assignment.filtered[vertex.order]

  # Identify which bi-clusters actually have vertices after disjoint assignment
  bc.present <- sort(unique(vertex.assignment.ordered))
  n.bc.present <- length(bc.present)

  cat("Bi-clusters with vertices after assignment:",
      paste(bc.present, collapse = ", "), "\n\n")

  # Create row annotation - only for bi-clusters that are present
  row.colors <- setNames(
    rainbow(n.bc)[bc.present],
    as.character(bc.present)
  )

  row.ha <- ComplexHeatmap::rowAnnotation(
    BiCluster = factor(vertex.assignment.ordered, levels = bc.present),
    col = list(BiCluster = row.colors),
    show_legend = !show.bicluster.labels,
    annotation_name_side = "top",
    simple_anno_size = grid::unit(5, "mm")
  )

  # Create column annotation for phylotypes (optional)
  if (show.phylotype.biclusters) {
    phylo.membership <- matrix(0, nrow = length(phylotype.indices), ncol = n.bc)
    for (i in 1:n.bc) {
      phylo.in.bc <- biclust.result@NumberxCol[i, phylotype.indices]
      phylo.membership[, i] <- phylo.in.bc * 1
    }
    colnames(phylo.membership) <- paste0("BC", 1:n.bc)

    col.ha <- ComplexHeatmap::HeatmapAnnotation(
      BiCluster = phylo.membership,
      col = list(BiCluster = circlize::colorRamp2(c(0, 1), c("white", "darkblue"))),
      show_legend = FALSE,
      annotation_name_side = "left",
      annotation_name_rot = 0
    )
  } else {
    col.ha <- NULL
  }

  # Color scheme for co-monotonicity
  col.fun <- circlize::colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    c("blue", "lightblue", "white", "pink", "red")
  )

  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    filtered.matrix,
    name = "Co-monotonicity",
    col = col.fun,
    cluster_rows = cluster.rows,
    cluster_columns = cluster.cols,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = "Vertices",
    row_title_side = "left",
    column_title = "Phylotypes",
    top_annotation = col.ha,
    right_annotation = row.ha,
    row_split = factor(vertex.assignment.ordered, levels = bc.present),
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    row_gap = grid::unit(2, "mm"),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Co-mono",
      direction = "horizontal"
    ),
    use_raster = nrow(filtered.matrix) > 1000 || ncol(filtered.matrix) > 1000,
    raster_quality = 2
  )

  # Add bi-cluster labels if requested
  # CRITICAL FIX: Only create labels for bi-clusters that are actually present
  if (show.bicluster.labels) {
    ht <- ht + ComplexHeatmap::rowAnnotation(
      BC = ComplexHeatmap::anno_block(
        gp = grid::gpar(fill = NA, col = NA),
        labels = bc.present,  # Use only present bi-clusters
        labels_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_rot = 0,
        which = "row"
      ),
      width = grid::unit(8, "mm")
    )
  }

  ht.drawn <- ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")

  return(invisible(list(
    heatmap = ht,
    drawn = ht.drawn,
    filtered.matrix = filtered.matrix,
    vertex.assignment = vertex.assignment,
    phylotype.indices = phylotype.indices,
    vertex.indices = vertex.indices[vertex.order],
    n.vertices.removed = n.vertices.removed,
    n.phylotypes.removed = n.phylotypes.removed,
    biclusters.present = bc.present
  )))
}

#' Summarize BCs4vd Bi-clustering Results
#'
#' Generates a summary table of bi-clusters identified by the BCs4vd sparse
#' singular value decomposition method, providing statistics on the size and
#' composition of each bi-cluster in terms of vertices (graph samples) and
#' phylotypes (bacterial taxa).
#'
#' @param biclust.result An object of class \code{Biclust} containing bi-clustering
#'   results from the BCs4vd method. Must contain slots \code{@Number} (number of
#'   bi-clusters), \code{@RowxNumber} (n.vertices x n.biclusters logical matrix
#'   indicating vertex membership), and \code{@NumberxCol} (n.biclusters x n.phylotypes
#'   logical matrix indicating phylotype membership).
#' @param base.adj.list Optional list of adjacency lists representing the base graph
#'   structure. If provided, the function will compute the number of connected components
#'   for the vertex set in each bi-cluster. Each element should be an integer vector
#'   of neighbor indices (1-based). Default is \code{NULL}.
#' @param vertex.names Optional character vector of vertex names or identifiers.
#'   Length must equal the number of rows in the original data matrix. Default is
#'   \code{NULL}.
#' @param phylotype.names Optional character vector of phylotype names. Length must
#'   equal the number of columns in the original data matrix. If \code{NULL}, phylotype
#'   indices will be used. Default is \code{NULL}.
#'
#' @return A data frame with one row per bi-cluster containing the following columns:
#'   \item{BiCluster}{Integer bi-cluster identifier (1 to n.biclusters)}
#'   \item{N_Vertices}{Number of vertices (graph samples) in the bi-cluster}
#'   \item{N_Phylotypes}{Number of phylotypes (bacterial taxa) in the bi-cluster}
#'   \item{Block_Size}{Total number of matrix cells in the rectangular block
#'     (N_Vertices x N_Phylotypes)}
#'   \item{N_Components}{Number of connected components in the base graph for this
#'     vertex set (only if \code{base.adj.list} is provided)}
#'   \item{Disconnected}{Logical indicating whether the vertex set spans multiple
#'     disconnected components (only if \code{base.adj.list} is provided)}
#'   \item{Vertex_Indices}{Comma-separated string of vertex indices (if requested)}
#'   \item{Phylotype_Indices}{Comma-separated string of phylotype indices (if requested)}
#'
#' @details
#' This function provides a comprehensive summary of bi-clustering results for
#' co-monotonicity matrices where rows represent vertices in a graph (samples or
#' spatial locations) and columns represent phylotypes (bacterial taxa). Each
#' bi-cluster identifies a rectangular block in the matrix where a subset of
#' phylotypes exhibits coordinated co-monotonicity patterns with the response
#' variable across a subset of vertices.
#'
#' When a base graph adjacency list is provided, the function checks whether the
#' vertices in each bi-cluster form a connected subgraph. Disconnected bi-clusters
#' are particularly interesting as they indicate that functionally similar regions
#' (similar co-monotonicity patterns) are geometrically separated in the sample
#' space. This can reveal convergent ecological patterns or redundant functional
#' pathways operating in different parts of the compositional manifold.
#'
#' The block size (N_Vertices x N_Phylotypes) provides a measure of the overall
#' size of each bi-cluster and can be used to prioritize bi-clusters for further
#' analysis. Larger blocks may represent more robust patterns, while smaller blocks
#' might capture more specific or localized associations.
#'
#' @note This function is specifically designed for bi-clustering results from the
#'   BCs4vd method applied to co-monotonicity matrices in the context of geometric
#'   data analysis for microbiome research. The interpretation assumes that rows
#'   represent samples arranged on a graph and columns represent features (phylotypes).
#'
#' @seealso \code{\link[s4vd]{BCs4vd}}, \code{\link[biclust]{biclust}}
#'
#' @examples
#' \dontrun{
#' library(s4vd)
#'
#' # Perform bi-clustering
#' bc.result <- biclust(comono.smooth,
#'                      method = BCs4vd(),
#'                      steps = 100,
#'                      pceru = 0.05,
#'                      pcerv = 0.05,
#'                      nbiclust = 15)
#'
#' # Basic summary without connectivity checking
#' bc.summary <- summarize.bcs4vd.biclusters(bc.result)
#' print(bc.summary)
#'
#' # Summary with connectivity analysis
#' bc.summary <- summarize.bcs4vd.biclusters(bc.result,
#'                                           base.adj.list = graph.adj.list,
#'                                           phylotype.names = colnames(comono.smooth))
#'
#' # Find disconnected bi-clusters
#' disconnected <- bc.summary[bc.summary$Disconnected, ]
#' cat("Found", nrow(disconnected), "disconnected bi-clusters\n")
#'
#' # Find large bi-clusters
#' large <- bc.summary[bc.summary$Block_Size > 1000, ]
#' }
#'
#' @export
summarize.bcs4vd.biclusters <- function(biclust.result,
                                        base.adj.list = NULL,
                                        vertex.names = NULL,
                                        phylotype.names = NULL) {

  n.bc <- biclust.result@Number

  if (n.bc == 0) {
    warning("No bi-clusters found in biclust.result")
    return(data.frame())
  }

  # Initialize summary data frame
  summary.df <- data.frame(
    BiCluster = 1:n.bc,
    N_Vertices = integer(n.bc),
    N_Phylotypes = integer(n.bc),
    Block_Size = integer(n.bc)
  )

  # Compute basic statistics
  for (i in 1:n.bc) {
    summary.df$N_Vertices[i] <- sum(biclust.result@RowxNumber[, i])
    summary.df$N_Phylotypes[i] <- sum(biclust.result@NumberxCol[i, ])
    summary.df$Block_Size[i] <- summary.df$N_Vertices[i] * summary.df$N_Phylotypes[i]
  }

  # Add connectivity information if base graph provided
  if (!is.null(base.adj.list)) {
    summary.df$N_Components <- integer(n.bc)
    summary.df$Disconnected <- logical(n.bc)

    for (i in 1:n.bc) {
      vertices <- which(biclust.result@RowxNumber[, i])
      n.comp <- check_connected_components(base.adj.list, vertices)
      summary.df$N_Components[i] <- n.comp
      summary.df$Disconnected[i] <- n.comp > 1
    }
  }

  # Add vertex indices as comma-separated string (optional, for small sets)
  summary.df$Vertex_Indices <- sapply(1:n.bc, function(i) {
    vertices <- which(biclust.result@RowxNumber[, i])
    if (length(vertices) <= 20) {
      paste(vertices, collapse = ", ")
    } else {
      paste0(paste(vertices[1:10], collapse = ", "), ", ..., ",
             paste(vertices[(length(vertices)-4):length(vertices)], collapse = ", "))
    }
  })

  # Add phylotype indices/names
  summary.df$Phylotype_Indices <- sapply(1:n.bc, function(i) {
    phylotypes <- which(biclust.result@NumberxCol[i, ])
    if (!is.null(phylotype.names)) {
      names <- phylotype.names[phylotypes]
      if (length(names) <= 20) {
        paste(names, collapse = ", ")
      } else {
        paste0(paste(names[1:10], collapse = ", "), ", ..., ",
               paste(names[(length(names)-4):length(names)], collapse = ", "))
      }
    } else {
      if (length(phylotypes) <= 20) {
        paste(phylotypes, collapse = ", ")
      } else {
        paste0(paste(phylotypes[1:10], collapse = ", "), ", ..., ",
               paste(phylotypes[(length(phylotypes)-4):length(phylotypes)], collapse = ", "))
      }
    }
  })

  # Print summary
  cat("Bi-cluster Summary Statistics\n")
  cat("=============================\n")
  cat("Total bi-clusters:", n.bc, "\n")
  cat("Vertices per bi-cluster: min =", min(summary.df$N_Vertices),
      ", median =", median(summary.df$N_Vertices),
      ", max =", max(summary.df$N_Vertices), "\n")
  cat("Phylotypes per bi-cluster: min =", min(summary.df$N_Phylotypes),
      ", median =", median(summary.df$N_Phylotypes),
      ", max =", max(summary.df$N_Phylotypes), "\n")
  cat("Block sizes: min =", min(summary.df$Block_Size),
      ", median =", median(summary.df$Block_Size),
      ", max =", max(summary.df$Block_Size), "\n")

  if (!is.null(base.adj.list)) {
    n.disconnected <- sum(summary.df$Disconnected)
    cat("\nConnectivity Analysis:\n")
    cat("  Disconnected bi-clusters:", n.disconnected, "(",
        round(100 * n.disconnected / n.bc, 1), "%)\n")
    cat("  Connected bi-clusters:", n.bc - n.disconnected, "(",
        round(100 * (n.bc - n.disconnected) / n.bc, 1), "%)\n")

    if (n.disconnected > 0) {
      cat("  Max components in a bi-cluster:",
          max(summary.df$N_Components[summary.df$Disconnected]), "\n")
    }
  }

  cat("\n")

  return(summary.df)
}

#' Visualize a Single BCs4vd Bi-cluster
#'
#' Creates a detailed heatmap visualization of a specific bi-cluster identified by
#' the BCs4vd sparse singular value decomposition method, showing the rectangular
#' block of co-monotonicity coefficients for the selected vertices and phylotypes.
#' Optionally analyzes the connectivity structure of the vertex set in the base graph
#' to determine whether the bi-cluster represents a geometrically coherent region.
#'
#' @param comono.matrix Numeric matrix of co-monotonicity coefficients with dimensions
#'   n.vertices x n.phylotypes. Rows represent vertices in the graph (samples or
#'   spatial locations) and columns represent phylotypes (bacterial taxa). Values
#'   should typically be in the range \code{[-1, 1]}, where positive values indicate
#'   similar directional associations with the response variable and negative values
#'   indicate opposite associations.
#' @param biclust.result An object of class \code{Biclust} containing bi-clustering
#'   results from the BCs4vd method. Must contain slots \code{@Number} (total number
#'   of bi-clusters), \code{@RowxNumber} (vertex membership matrix), and
#'   \code{@NumberxCol} (phylotype membership matrix).
#' @param bc.number Integer specifying which bi-cluster to visualize (must be between
#'   1 and \code{biclust.result@Number}).
#' @param base.adj.list Optional list of adjacency lists representing the base graph
#'   structure. If provided, the function performs connectivity analysis to determine
#'   whether the vertices in this bi-cluster form a connected subgraph. Each element
#'   should be an integer vector of neighbor indices (1-based). Default is \code{NULL}.
#' @param show.row.names Logical indicating whether to display row (vertex) names on
#'   the heatmap. Only recommended for small bi-clusters. Default is \code{FALSE}.
#' @param show.col.names Logical indicating whether to display column (phylotype)
#'   names on the heatmap. Default is \code{TRUE}.
#'
#' @return Returns invisibly a list containing:
#'   \item{heatmap}{The \code{Heatmap} object from ComplexHeatmap package}
#'   \item{drawn}{The drawn heatmap object (output of \code{draw()})}
#'   \item{vertices}{Integer vector of vertex indices in this bi-cluster}
#'   \item{phylotypes}{Integer vector of phylotype indices in this bi-cluster}
#'   \item{n.components}{Number of connected components in the base graph for this
#'     vertex set (only if \code{base.adj.list} provided)}
#'   \item{disconnected}{Logical indicating whether the vertex set spans multiple
#'     disconnected components (only if \code{base.adj.list} provided)}
#'   \item{matrix}{The extracted rectangular block matrix}
#'
#' @details
#' This function extracts and visualizes the rectangular block corresponding to a
#' single bi-cluster from a co-monotonicity matrix. The extracted block has
#' dimensions n.vertices.in.cluster x n.phylotypes.in.cluster and shows the
#' co-monotonicity patterns between the selected phylotypes and the response
#' variable across the selected vertices.
#'
#' The heatmap uses hierarchical clustering to reorder both rows and columns,
#' revealing internal structure within the bi-cluster. A diverging color scheme
#' (blue-white-red) represents negative, neutral, and positive co-monotonicity
#' values respectively.
#'
#' When a base graph adjacency list is provided, the function performs breadth-first
#' search to determine the number of connected components among the vertices in the
#' bi-cluster. A disconnected bi-cluster (multiple components) indicates that
#' phylotypes exhibit similar co-monotonicity patterns in geometrically separated
#' regions of the sample space. This can reveal important biological insights about
#' convergent ecological processes or redundant functional pathways operating
#' independently in different parts of the compositional manifold.
#'
#' The function prints diagnostic information to the console including the number
#' of vertices, phylotypes, and connectivity status (if applicable).
#'
#' @note This function requires the ComplexHeatmap and circlize packages. It is
#'   specifically designed for bi-clustering results from sparse SVD methods applied
#'   to co-monotonicity matrices in gradient flow analysis. For visualization of
#'   all bi-clusters simultaneously, see \code{\link{visualize.bcs4vd.biclusters}}.
#'
#' @seealso \code{\link{visualize.bcs4vd.biclusters}}, \code{\link{summarize.bcs4vd.biclusters}},
#'   \code{\link[s4vd]{BCs4vd}}
#'
#' @examples
#' \dontrun{
#' library(s4vd)
#' library(ComplexHeatmap)
#'
#' # Perform bi-clustering
#' bc.result <- biclust(comono.smooth,
#'                      method = BCs4vd(),
#'                      steps = 100,
#'                      pceru = 0.05,
#'                      pcerv = 0.05,
#'                      nbiclust = 15)
#'
#' # Visualize first bi-cluster without connectivity analysis
#' bc1 <- visualize.bcs4vd.bicluster(comono.smooth, bc.result, bc.number = 1)
#'
#' # Visualize with connectivity analysis
#' bc1 <- visualize.bcs4vd.bicluster(comono.smooth, bc.result,
#'                                   bc.number = 1,
#'                                   base.adj.list = graph.adj.list)
#'
#' # Access the extracted data
#' cat("Vertices:", length(bc1$vertices), "\n")
#' cat("Phylotypes:", length(bc1$phylotypes), "\n")
#' cat("Disconnected:", bc1$disconnected, "\n")
#'
#' # Visualize multiple bi-clusters
#' for (i in 1:5) {
#'   visualize.bcs4vd.bicluster(comono.smooth, bc.result, bc.number = i)
#' }
#' }
#'
#' @export
visualize.bcs4vd.bicluster <- function(comono.matrix,
                                       biclust.result,
                                       bc.number,
                                       base.adj.list = NULL,
                                       show.row.names = FALSE,
                                       show.col.names = TRUE) {

  # Validate inputs
  if (bc.number < 1 || bc.number > biclust.result@Number) {
    stop(paste("bc.number must be between 1 and", biclust.result@Number))
  }

  if (nrow(comono.matrix) != nrow(biclust.result@RowxNumber)) {
    stop("Number of rows in comono.matrix does not match biclust.result")
  }

  if (ncol(comono.matrix) != ncol(biclust.result@NumberxCol)) {
    stop("Number of columns in comono.matrix does not match biclust.result")
  }

  # Extract vertex indices (rows) and phylotype indices (columns)
  vertices <- which(biclust.result@RowxNumber[, bc.number])
  phylotypes <- which(biclust.result@NumberxCol[bc.number, ])

  cat("Bi-cluster", bc.number, ":\n")
  cat("  Vertices (samples):", length(vertices), "\n")
  cat("  Phylotypes (taxa):", length(phylotypes), "\n")

  # Initialize return list
  result.list <- list(
    vertices = vertices,
    phylotypes = phylotypes
  )

  # Check connectivity if base graph provided
  if (!is.null(base.adj.list)) {
    n.comp <- check_connected_components(base.adj.list, vertices)
    result.list$n.components <- n.comp
    result.list$disconnected <- n.comp > 1

    cat("  Connected components in base graph:", n.comp, "\n")
    if (n.comp > 1) {
      cat("  -> DISCONNECTED: Functionally similar but geometrically separated!\n")
    } else {
      cat("  -> CONNECTED: Functionally and geometrically coherent\n")
    }
  }
  cat("\n")

  # Extract the rectangular block
  block <- comono.matrix[vertices, phylotypes, drop = FALSE]
  result.list$matrix <- block

  # Create phylotype labels
  if (show.col.names && !is.null(colnames(comono.matrix))) {
    phylo.labels <- colnames(comono.matrix)[phylotypes]
  } else if (show.col.names) {
    phylo.labels <- paste0("P", phylotypes)
  } else {
    phylo.labels <- NULL
  }

  # Create vertex labels
  if (show.row.names && !is.null(rownames(comono.matrix))) {
    vertex.labels <- rownames(comono.matrix)[vertices]
  } else if (show.row.names) {
    vertex.labels <- paste0("V", vertices)
  } else {
    vertex.labels <- NULL
  }

  # Color scheme for co-monotonicity
  col.fun <- circlize::colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    c("blue", "lightblue", "white", "pink", "red")
  )

  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    block,
    name = "Co-mono",
    col = col.fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = show.row.names,
    show_column_names = show.col.names,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    column_labels = phylo.labels,
    row_labels = vertex.labels,
    column_title = paste0("Bi-cluster ", bc.number, ": ",
                         length(phylotypes), " phylotypes"),
    row_title = paste0(length(vertices), " vertices"),
    column_names_rot = 45,
    use_raster = length(vertices) > 500,
    raster_quality = 2,
    heatmap_legend_param = list(
      title = "Co-mono",
      direction = "vertical"
    )
  )

  # Draw heatmap
  ht.drawn <- ComplexHeatmap::draw(ht, heatmap_legend_side = "right")

  result.list$heatmap <- ht
  result.list$drawn <- ht.drawn

  return(invisible(result.list))
}

# Helper function for connectivity checking (include if not already defined elsewhere)
#' Check Connected Components in Subgraph
#'
#' Internal helper function to determine the number of connected components
#' in a subgraph induced by a subset of vertices.
#'
#' @param adj.list List of adjacency lists for the full graph (1-based indexing)
#' @param vertex.indices Integer vector of vertex indices to include in subgraph
#'
#' @return Integer number of connected components
#' @keywords internal
check_connected_components <- function(adj.list, vertex.indices) {
  n <- length(vertex.indices)

  if (n == 0) return(0)
  if (n == 1) return(1)

  # Build index mapping
  idx.map <- setNames(1:n, vertex.indices)

  # BFS to find components
  visited <- rep(FALSE, n)
  n.components <- 0

  for (start in 1:n) {
    if (!visited[start]) {
      n.components <- n.components + 1

      # BFS from this vertex
      queue <- start
      visited[start] <- TRUE

      while (length(queue) > 0) {
        current <- queue[1]
        queue <- queue[-1]

        # Get neighbors in original graph
        orig.idx <- vertex.indices[current]

        # Check if this vertex exists in adj.list
        if (orig.idx <= length(adj.list) && !is.null(adj.list[[orig.idx]])) {
          neighbors.orig <- adj.list[[orig.idx]]

          # Find neighbors in our vertex set
          neighbors.in.set <- neighbors.orig[neighbors.orig %in% vertex.indices]

          if (length(neighbors.in.set) > 0) {
            neighbors.local <- idx.map[as.character(neighbors.in.set)]

            for (neighbor in neighbors.local) {
              if (!visited[neighbor]) {
                visited[neighbor] <- TRUE
                queue <- c(queue, neighbor)
              }
            }
          }
        }
      }
    }
  }

  return(n.components)
}

#' Compute Phylotype Frequency Across BCs4vd Bi-clusters
#'
#' Generates a frequency table showing how many bi-clusters each phylotype appears
#' in, sorted in decreasing order of frequency. This analysis identifies hub
#' phylotypes that participate in multiple co-monotonicity patterns across different
#' regions of the sample space, as well as specific phylotypes that exhibit
#' localized associations with the response variable.
#'
#' @param biclust.result An object of class \code{Biclust} containing bi-clustering
#'   results from the BCs4vd method. Must contain slots \code{@Number} (number of
#'   bi-clusters) and \code{@NumberxCol} (n.biclusters x n.phylotypes logical matrix
#'   indicating phylotype membership in each bi-cluster).
#' @param phylo.names Optional character vector of phylotype names. Length must equal
#'   the number of phylotypes (columns) in the original co-monotonicity matrix. If
#'   \code{NULL}, phylotype indices will be used as identifiers. Default is \code{NULL}.
#'
#' @return A data frame with one row per phylotype that appears in at least one
#'   bi-cluster, containing the following columns:
#'   \item{Phylotype_Index}{Integer index of the phylotype in the original matrix}
#'   \item{N_BiClusters}{Number of bi-clusters containing this phylotype}
#'   \item{Phylotype_Name}{Name of the phylotype (only if \code{phylo.names} provided)}
#'
#'   The data frame is sorted by \code{N_BiClusters} in decreasing order. Returns
#'   \code{NULL} if no phylotypes are found in any bi-clusters.
#'
#' @details
#' This function provides a summary of phylotype participation across bi-clusters,
#' which is useful for identifying different classes of phylotypes based on their
#' co-monotonicity patterns:
#'
#' Hub phylotypes appear in many bi-clusters, indicating that their associations
#' with the response variable (such as spontaneous preterm birth risk) are consistent
#' across multiple regions of the sample space. These phylotypes may represent core
#' members of the microbial community whose functional roles are broadly conserved
#' across different ecological contexts or patient subpopulations.
#'
#' Specific phylotypes appear in few bi-clusters (ideally just one), indicating
#' localized or context-dependent associations with the response variable. These
#' phylotypes may be sensitive to specific environmental conditions or may interact
#' with other community members in region-specific ways. They represent potential
#' biomarkers for distinguishing between different ecological states or patient
#' subgroups.
#'
#' Phylotypes that appear in exactly zero bi-clusters are excluded from the output
#' table. These phylotypes may have weak or inconsistent co-monotonicity patterns
#' that did not meet the sparsity criteria during bi-clustering, or they may
#' exhibit high variability that prevents stable cluster assignment.
#'
#' The function prints diagnostic information to the console including the total
#' number of phylotypes, how many appear in bi-clusters, and a frequency distribution
#' showing how many phylotypes appear in 1, 2, 3, etc. bi-clusters.
#'
#' @note This function is designed specifically for BCs4vd bi-clustering results
#'   applied to co-monotonicity matrices in microbiome research. For a more detailed
#'   analysis including which specific bi-clusters each phylotype belongs to, see
#'   \code{phylotype.bcs4vd.bicluster.detailed}.
#'
#' @seealso \code{phylotype.bcs4vd.bicluster.detailed},
#'   \code{\link{summarize.bcs4vd.biclusters}}, \code{\link[s4vd]{BCs4vd}}
#'
#' @examples
#' \dontrun{
#' library(s4vd)
#'
#' # Perform bi-clustering
#' bc.result <- biclust(comono.smooth,
#'                      method = BCs4vd(),
#'                      steps = 100,
#'                      pceru = 0.05,
#'                      pcerv = 0.05,
#'                      nbiclust = 15)
#'
#' # Get phylotype frequency table with indices only
#' phylo.freq <- phylotype.bcs4vd.bicluster.frequency(bc.result)
#' head(phylo.freq)
#'
#' # Get phylotype frequency table with names
#' phylo.freq <- phylotype.bcs4vd.bicluster.frequency(
#'   bc.result,
#'   phylo.names = colnames(comono.smooth)
#' )
#'
#' # Identify hub phylotypes (appear in many bi-clusters)
#' hub.phylotypes <- phylo.freq[phylo.freq$N_BiClusters >= 5, ]
#' cat("Hub phylotypes:\n")
#' print(hub.phylotypes)
#'
#' # Identify specific phylotypes (appear in only one bi-cluster)
#' specific.phylotypes <- phylo.freq[phylo.freq$N_BiClusters == 1, ]
#' cat("Specific phylotypes:", nrow(specific.phylotypes), "\n")
#'
#' # Plot frequency distribution
#' barplot(table(phylo.freq$N_BiClusters),
#'         xlab = "Number of Bi-clusters",
#'         ylab = "Number of Phylotypes",
#'         main = "Phylotype Bi-cluster Frequency Distribution")
#' }
#'
#' @export
phylotype.bcs4vd.bicluster.frequency <- function(biclust.result, phylo.names = NULL) {

  # NumberxCol is n.biclusters x n.phylotypes
  # Each column represents a phylotype
  # Each row represents a bi-cluster
  # Value is TRUE/1 if phylotype is in that bi-cluster

  n.phylotypes <- ncol(biclust.result@NumberxCol)
  n.biclusters <- nrow(biclust.result@NumberxCol)

  # Count how many bi-clusters each phylotype appears in
  phylo.counts <- colSums(biclust.result@NumberxCol)

  # Keep only phylotypes that appear in at least one bi-cluster
  phylotypes.in.clusters <- which(phylo.counts > 0)

  if (length(phylotypes.in.clusters) == 0) {
    cat("No phylotypes found in any bi-clusters\n")
    return(NULL)
  }

  # Create frequency table
  freq.table <- data.frame(
    Phylotype_Index = phylotypes.in.clusters,
    N_BiClusters = phylo.counts[phylotypes.in.clusters]
  )

  # Add phylotype names if provided
  if (!is.null(phylo.names)) {
    if (length(phylo.names) != n.phylotypes) {
      warning("Length of phylo.names does not match number of phylotypes. Names not added.")
    } else {
      freq.table$Phylotype_Name <- phylo.names[phylotypes.in.clusters]
    }
  }

  # Sort by frequency (decreasing)
  freq.table <- freq.table[order(freq.table$N_BiClusters, decreasing = TRUE), ]

  # Reset row names
  rownames(freq.table) <- NULL

  # Print summary
  cat("Phylotype Bi-cluster Frequency Summary\n")
  cat("======================================\n")
  cat("Total phylotypes in dataset:", n.phylotypes, "\n")
  cat("Phylotypes in at least 1 bi-cluster:", length(phylotypes.in.clusters),
      "(", round(100 * length(phylotypes.in.clusters) / n.phylotypes, 1), "%)\n")
  cat("Phylotypes in multiple bi-clusters:", sum(phylo.counts > 1),
      "(", round(100 * sum(phylo.counts > 1) / n.phylotypes, 1), "%)\n")
  cat("Max bi-clusters per phylotype:", max(phylo.counts), "\n")
  cat("\nFrequency distribution:\n")
  freq.dist <- table(phylo.counts[phylo.counts > 0])
  print(freq.dist)
  cat("\n")

  return(freq.table)
}

#' Diagnose Bi-cluster Coherence
#'
#' Checks if a bi-cluster satisfies the rank-one assumption by examining
#' pattern heterogeneity within the block.
#'
#' @param data.matrix Original data matrix
#' @param biclust.result Biclust object
#' @param bc.number Which bi-cluster to diagnose
#'
#' @return List with coherence diagnostics
#'
diagnose.bicluster.coherence <- function(data.matrix, biclust.result, bc.number) {

  # Extract bi-cluster
  rows <- which(biclust.result@RowxNumber[, bc.number])
  cols <- which(biclust.result@NumberxCol[bc.number, ])

  block <- data.matrix[rows, cols, drop = FALSE]

  # SVD of the block
  svd.block <- svd(block)

  # Variance explained by first singular value (should be high for rank-1)
  var.explained <- svd.block$d^2 / sum(svd.block$d^2)

  # Within-block correlation structure
  # If rank-1, all row vectors should be highly correlated (up to sign)
  row.cors <- cor(t(block))
  row.cors.abs <- abs(row.cors[upper.tri(row.cors)])

  col.cors <- cor(block)
  col.cors.abs <- abs(col.cors[upper.tri(col.cors)])

  # Sign consistency - check if values are consistently positive or negative
  # For each phylotype, what proportion of vertices have same sign?
  sign.consistency <- apply(block, 2, function(col) {
    if (all(col == 0)) return(NA)
    pos <- sum(col > 0)
    neg <- sum(col < 0)
    max(pos, neg) / length(col)
  })

  # Hierarchical clustering of rows to detect subgroups
  if (nrow(block) > 2) {
    row.dist <- dist(block)
    row.hclust <- hclust(row.dist, method = "ward.D2")
    # Try 2-cluster split
    row.clusters <- cutree(row.hclust, k = 2)

    # Silhouette to see if split is meaningful
    sil <- cluster::silhouette(row.clusters, row.dist)
    avg.sil <- mean(sil[, 3])
  } else {
    row.clusters <- NULL
    avg.sil <- NA
  }

  cat("Bi-cluster", bc.number, "Coherence Diagnostics:\n")
  cat("=====================================\n")
  cat("Block size:", nrow(block), "x", ncol(block), "\n")
  cat("Variance explained by rank-1:", round(var.explained[1] * 100, 1), "%\n")
  cat("First 5 singular values:", round(svd.block$d[1:min(5, length(svd.block$d))], 2), "\n")
  cat("\nCorrelation structure:\n")
  cat("  Median abs row correlation:", round(median(row.cors.abs, na.rm = TRUE), 3), "\n")
  cat("  Median abs col correlation:", round(median(col.cors.abs, na.rm = TRUE), 3), "\n")
  cat("\nSign consistency (per phylotype):\n")
  cat("  Mean:", round(mean(sign.consistency, na.rm = TRUE), 3), "\n")
  cat("  Median:", round(median(sign.consistency, na.rm = TRUE), 3), "\n")
  cat("  Min:", round(min(sign.consistency, na.rm = TRUE), 3), "\n")

  if (!is.na(avg.sil)) {
    cat("\nSubstructure detection:\n")
    cat("  2-cluster silhouette:", round(avg.sil, 3), "\n")
    if (avg.sil > 0.25) {
      cat("  -> Significant substructure detected! Consider splitting.\n")
    }
  }
  cat("\n")

  return(invisible(list(
    block = block,
    var.explained = var.explained,
    row.correlations = row.cors,
    col.correlations = col.cors,
    sign.consistency = sign.consistency,
    row.clusters = row.clusters,
    avg.silhouette = avg.sil,
    svd = svd.block
  )))
}

#' Recursive Bi-clustering Refinement with Enhanced Diagnostics
#'
#' Recursively applies bi-clustering to refine heterogeneous bi-clusters.
#' Includes multiple refinement strategies and detailed convergence diagnostics.
#'
#' @param data.matrix Original data matrix
#' @param biclust.result Initial bi-clustering result
#' @param bc.to.refine Which bi-clusters to refine (indices)
#' @param coherence.threshold Variance explained threshold for refinement
#' @param pceru Per-comparison error rate for rows
#' @param pcerv Per-comparison error rate for columns
#' @param max.steps Maximum BCs4vd iterations. Default 500 (increased from 100)
#' @param nbiclust Maximum number of sub-clusters to find
#' @param fallback.method Method to use if BCs4vd fails: "kmeans", "hclust", or "none"
#' @param fallback.k Number of clusters for fallback method
#' @param verbose Print detailed diagnostics
#'
#' @return List of refined bi-clusters with metadata
#'
refine.biclusters <- function(data.matrix,
                              biclust.result,
                              bc.to.refine = NULL,
                              coherence.threshold = 0.7,
                              pceru = 0.1,
                              pcerv = 0.1,
                              max.steps = 500,
                              nbiclust = 3,
                              fallback.method = c("kmeans", "hclust", "none"),
                              fallback.k = 2,
                              verbose = TRUE) {

  fallback.method <- match.arg(fallback.method)
  pkg.biclust <- "biclust"
  pkg.s4vd <- "s4vd"
  if (!requireNamespace(pkg.biclust, quietly = TRUE) ||
      !requireNamespace(pkg.s4vd, quietly = TRUE)) {
    stop("refine.biclusters() requires optional packages 'biclust' and 's4vd'")
  }
  biclust_fn <- getExportedValue(pkg.biclust, "biclust")
  bcs4vd_fn <- getExportedValue(pkg.s4vd, "BCs4vd")

  # If not specified, identify bi-clusters needing refinement
  if (is.null(bc.to.refine)) {
    bc.to.refine <- c()
    for (i in 1:biclust.result@Number) {
      rows <- which(biclust.result@RowxNumber[, i])
      cols <- which(biclust.result@NumberxCol[i, ])  # Fixed: was bc.number
      block <- data.matrix[rows, cols, drop = FALSE]

      svd.block <- svd(block)
      var.exp <- (svd.block$d[1]^2) / sum(svd.block$d^2)

      if (var.exp < coherence.threshold && nrow(block) >= 8) {
        bc.to.refine <- c(bc.to.refine, i)
      }
    }
  }

  if (verbose) {
    cat("Refining bi-clusters:", paste(bc.to.refine, collapse = ", "), "\n")
    cat("PCER: rows =", pceru, ", cols =", pcerv, "\n")
    cat("Max steps:", max.steps, "\n\n")
  }

  refined.list <- list()

  for (bc.id in bc.to.refine) {
    if (verbose) {
      cat("Refining bi-cluster", bc.id, "...\n")
    }

    # Extract vertices and phylotypes
    rows <- which(biclust.result@RowxNumber[, bc.id])
    cols <- which(biclust.result@NumberxCol[bc.id, ])

    # Subset data
    sub.data <- data.matrix[rows, cols, drop = FALSE]

    if (verbose) {
      cat("  Sub-block size:", nrow(sub.data), "x", ncol(sub.data), "\n")
    }

    # Check if block is too small
    if (nrow(sub.data) < 8 || ncol(sub.data) < 4) {
      if (verbose) {
        cat("  WARNING: Block too small for reliable refinement\n\n")
      }
      next
    }

    # Apply bi-clustering to subset with increased steps
    sub.bc <- tryCatch({
      biclust_fn(sub.data,
             method = bcs4vd_fn(),
             steps = max.steps,
             pceru = pceru,
             pcerv = pcerv,
             nbiclust = nbiclust,
             row.min = 4,
             col.min = 4)
    }, error = function(e) {
      if (verbose) {
        cat("  ERROR in BCs4vd:", e$message, "\n")
      }
      return(NULL)
    })

    # Check convergence
    converged <- !is.null(sub.bc) && sub.bc@Number > 0

    if (converged) {
      if (verbose) {
        cat("  SUCCESS: Found", sub.bc@Number, "sub-clusters\n")
      }

      # Store refined results
      for (j in 1:sub.bc@Number) {
        sub.rows <- rows[which(sub.bc@RowxNumber[, j])]
        sub.cols <- cols[which(sub.bc@NumberxCol[j, ])]

        refined.list[[length(refined.list) + 1]] <- list(
          parent = bc.id,
          sub.cluster = j,
          rows = sub.rows,
          cols = sub.cols,
          method = "BCs4vd",
          converged = TRUE
        )
      }
    } else {
      # BCs4vd failed - try fallback method
      if (verbose) {
        cat("  WARNING: BCs4vd did not converge\n")
      }

      if (fallback.method != "none") {
        if (verbose) {
          cat("  Attempting fallback method:", fallback.method, "\n")
        }

        fallback.result <- apply.fallback.split(
          sub.data, rows, cols,
          fallback.method, fallback.k,
          verbose
        )

        if (!is.null(fallback.result)) {
          refined.list <- c(refined.list, fallback.result)
        }
      }
    }

    if (verbose) cat("\n")
  }

  if (verbose) {
    cat("Refinement complete: generated", length(refined.list), "sub-clusters\n")
  }

  return(refined.list)
}

#' Apply Fallback Splitting Method
#'
#' Uses k-means or hierarchical clustering when BCs4vd fails
#'
#' @keywords internal
apply.fallback.split <- function(sub.data, rows, cols, method, k, verbose) {

  # Detect natural number of clusters if k not specified
  if (is.null(k) || k < 2) {
    # Use gap statistic or silhouette
    sil.scores <- numeric(min(5, nrow(sub.data) - 1))
    for (test.k in 2:min(5, nrow(sub.data) - 1)) {
      if (method == "kmeans") {
        km.test <- kmeans(sub.data, centers = test.k, nstart = 25)
        sil <- cluster::silhouette(km.test$cluster, dist(sub.data))
      } else {
        hc.test <- hclust(dist(sub.data), method = "ward.D2")
        clust.test <- cutree(hc.test, k = test.k)
        sil <- cluster::silhouette(clust.test, dist(sub.data))
      }
      sil.scores[test.k - 1] <- mean(sil[, 3])
    }
    k <- which.max(sil.scores) + 1
    if (verbose) {
      cat("    Auto-detected k =", k, "(silhouette =",
          round(max(sil.scores), 3), ")\n")
    }
  }

  # Apply chosen method
  result.list <- list()

  if (method == "kmeans") {
    set.seed(123)
    km <- kmeans(sub.data, centers = k, nstart = 50)
    clusters <- km$cluster

    if (verbose) {
      cat("    K-means split: sizes =", paste(table(clusters), collapse = ", "), "\n")
    }

  } else if (method == "hclust") {
    hc <- hclust(dist(sub.data), method = "ward.D2")
    clusters <- cutree(hc, k = k)

    if (verbose) {
      cat("    Hierarchical split: sizes =", paste(table(clusters), collapse = ", "), "\n")
    }
  }

  # Create refined entries for each cluster
  for (j in 1:k) {
    cluster.rows <- rows[clusters == j]

    # For columns, keep all that have reasonable signal in this row subset
    cluster.data <- sub.data[clusters == j, , drop = FALSE]
    col.means <- colMeans(abs(cluster.data))
    active.cols <- which(col.means > quantile(col.means, 0.25))
    cluster.cols <- cols[active.cols]

    if (length(cluster.rows) >= 4 && length(cluster.cols) >= 4) {
      result.list[[length(result.list) + 1]] <- list(
        parent = NA,  # Will be set by caller
        sub.cluster = j,
        rows = cluster.rows,
        cols = cluster.cols,
        method = method,
        converged = FALSE
      )
    }
  }

  # Update parent IDs
  parent.id <- result.list[[1]]$parent
  for (j in seq_along(result.list)) {
    result.list[[j]]$parent <- parent.id
  }

  return(result.list)
}

#' Diagnose Why Refinement Failed
#'
#' Analyzes a bi-cluster to understand why BCs4vd refinement didn't converge
#'
#' @param data.matrix Original data matrix
#' @param biclust.result Biclust object
#' @param bc.number Which bi-cluster to diagnose
#'
diagnose.refinement.failure <- function(data.matrix, biclust.result, bc.number) {

  rows <- which(biclust.result@RowxNumber[, bc.number])
  cols <- which(biclust.result@NumberxCol[bc.number, ])
  block <- data.matrix[rows, cols, drop = FALSE]

  cat("Refinement Failure Diagnostics for BC", bc.number, "\n")
  cat("==========================================\n\n")

  # Size check
  cat("Block dimensions:", nrow(block), "x", ncol(block), "\n")
  if (nrow(block) < 8) {
    cat("  WARNING: Too few rows for stable bi-clustering\n")
  }
  if (ncol(block) < 4) {
    cat("  WARNING: Too few columns for stable bi-clustering\n")
  }
  cat("\n")

  # Coherence check
  svd.block <- svd(block)
  var.exp <- svd.block$d^2 / sum(svd.block$d^2)
  cat("Variance structure:\n")
  cat("  Rank-1 explained:", round(var.exp[1] * 100, 1), "%\n")
  cat("  Rank-2 explained:", round(sum(var.exp[1:2]) * 100, 1), "%\n")
  cat("  Effective rank:", sum(var.exp > 0.01), "\n")

  if (var.exp[1] > 0.95) {
    cat("  NOTE: Very high rank-1 - block may be too coherent to split\n")
  }
  cat("\n")

  # Correlation structure
  row.cors <- cor(t(block))
  diag(row.cors) <- NA
  cat("Row correlation structure:\n")
  cat("  Mean:", round(mean(abs(row.cors), na.rm = TRUE), 3), "\n")
  cat("  Median:", round(median(abs(row.cors), na.rm = TRUE), 3), "\n")
  cat("  Range:", round(range(abs(row.cors), na.rm = TRUE), 3), "\n")

  if (median(abs(row.cors), na.rm = TRUE) > 0.9) {
    cat("  NOTE: Very high inter-row correlation - limited substructure\n")
  }
  cat("\n")

  # Detect potential splits via hierarchical clustering
  row.dist <- dist(block)
  row.hclust <- hclust(row.dist, method = "ward.D2")

  cat("Hierarchical clustering analysis:\n")
  for (k in 2:min(4, nrow(block) - 1)) {
    clust <- cutree(row.hclust, k = k)
    sil <- cluster::silhouette(clust, row.dist)
    avg.sil <- mean(sil[, 3])
    sizes <- table(clust)
    cat("  k =", k, ": silhouette =", round(avg.sil, 3),
        ", sizes =", paste(sizes, collapse = ", "))
    if (avg.sil > 0.5) {
      cat(" [GOOD]")
    } else if (avg.sil > 0.25) {
      cat(" [MODERATE]")
    } else {
      cat(" [WEAK]")
    }
    cat("\n")
  }
  cat("\n")

  # Recommendation
  cat("Recommendation:\n")
  if (var.exp[1] > 0.95 && median(abs(row.cors), na.rm = TRUE) > 0.9) {
    cat("  Block is highly coherent - may not need refinement\n")
  } else {
    cat("  Try: (1) Increase max.steps to 1000\n")
    cat("       (2) Adjust PCER: try pceru = 0.05 or 0.15\n")
    cat("       (3) Use fallback.method = 'kmeans' or 'hclust'\n")
  }
}

#' Refine Bi-cluster Using K-means
#'
#' Split heterogeneous bi-cluster using k-means on vertex patterns
#'
#' @param data.matrix Numeric matrix used for biclustering.
#' @param biclust.result Object containing the original biclustering result.
#' @param bc.number Index of the bicluster to refine.
#' @param k Number of k-means subclusters used for refinement.
#'
kmeans.refine.bicluster <- function(data.matrix, biclust.result, bc.number, k = 2) {

  rows <- which(biclust.result@RowxNumber[, bc.number])
  cols <- which(biclust.result@NumberxCol[bc.number, ])

  block <- data.matrix[rows, cols, drop = FALSE]

  # K-means on rows
  set.seed(123)
  km <- kmeans(block, centers = k, nstart = 25)

  # Visualize split
  row.ha <- ComplexHeatmap::rowAnnotation(
    KMeans = factor(km$cluster),
    col = list(KMeans = setNames(rainbow(k), as.character(1:k)))
  )

  ht <- ComplexHeatmap::Heatmap(
    block[order(km$cluster), ],
    name = "Co-mono",
    col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    row_split = factor(km$cluster[order(km$cluster)]),
    column_title = paste("BC", bc.number, "- K-means Split"),
    right_annotation = row.ha
  )

  ComplexHeatmap::draw(ht)

  return(list(
    clusters = km$cluster,
    rows.by.cluster = split(rows, km$cluster),
    block = block
  ))
}

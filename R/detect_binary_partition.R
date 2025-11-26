#' Detect Strong Binary Partition in Dendrogram
#'
#' Specialized method for detecting two very strong, well-separated clusters
#' characterized by large gap in merge heights and small within-cluster variance
#'
#' @param data.matrix Data matrix (vertices × features)
#' @param linkage.method Hierarchical clustering linkage
#' @param min.cluster.fraction Minimum fraction of data in each cluster
#' @param plot.diagnostics Whether to create diagnostic plots
#'
#' @return List with binary partition and diagnostics
#'
#' @details
#' This function specifically looks for the "strong two-cluster" pattern where:
#' 1. Two large clusters contain most of the data
#' 2. Within-cluster variance is much smaller than between-cluster variance
#' 3. The merge height for combining these clusters is much larger than
#'    the heights within each cluster
#'
#' @export
detect.strong.binary.partition <- function(data.matrix,
                                           linkage.method = "ward.D2",
                                           min.cluster.fraction = 0.1,
                                           plot.diagnostics = TRUE) {
  
  n <- nrow(data.matrix)
  data.dist <- dist(data.matrix, method = "euclidean")
  data.hclust <- hclust(data.dist, method = linkage.method)
  
  merge.heights <- data.hclust$height
  n.merges <- length(merge.heights)
  
  # For each possible k=2 cut, evaluate the partition quality
  # We'll try cutting at different heights and find the one with
  # maximum gap relative to within-cluster structure
  
  best.score <- -Inf
  best.partition <- NULL
  best.cut.height <- NULL
  
  cat("Searching for optimal binary partition...\n")
  
  # Try cutting at different percentiles of the height range
  height.range <- range(merge.heights)
  test.heights <- seq(
    quantile(merge.heights, 0.1),
    quantile(merge.heights, 0.9),
    length.out = 100
  )
  
  scores <- numeric(length(test.heights))
  
  for (i in seq_along(test.heights)) {
    cut.height <- test.heights[i]
    clusters <- cutree(data.hclust, h = cut.height)
    
    k <- length(unique(clusters))
    if (k != 2) next
    
    # Check size constraint
    cluster.sizes <- table(clusters)
    if (any(cluster.sizes < n * min.cluster.fraction)) next
    
    # Compute partition quality
    # 1. Within-cluster tightness
    within.var <- sum(sapply(1:2, function(cl) {
      if (sum(clusters == cl) <= 1) return(0)
      cl.data <- data.matrix[clusters == cl, ]
      cl.center <- colMeans(cl.data)
      sum(sweep(cl.data, 2, cl.center)^2)
    }))
    
    # 2. Between-cluster separation
    center1 <- colMeans(data.matrix[clusters == 1, ])
    center2 <- colMeans(data.matrix[clusters == 2, ])
    between.dist <- sqrt(sum((center1 - center2)^2))
    
    # Score: maximize separation relative to within-variance
    # This is essentially the F-statistic from ANOVA
    score <- between.dist^2 / (within.var / (n - 2))
    scores[i] <- score
    
    if (score > best.score) {
      best.score <- score
      best.partition <- clusters
      best.cut.height <- cut.height
    }
  }
  
  if (is.null(best.partition)) {
    stop("Could not find valid binary partition satisfying constraints")
  }
  
  cat("Found optimal binary partition at height:", 
      round(best.cut.height, 3), "\n")
  cat("Partition quality score:", round(best.score, 2), "\n\n")
  
  # Detailed diagnostics
  diagnostics <- compute.binary.partition.diagnostics(
    data.matrix, 
    best.partition, 
    data.hclust,
    best.cut.height
  )
  
  # Plotting
  if (plot.diagnostics) {
    plot.binary.partition.diagnostics(
      data.matrix,
      best.partition,
      data.hclust,
      best.cut.height,
      scores,
      test.heights,
      diagnostics
    )
  }
  
  return(list(
    clusters = best.partition,
    cut.height = best.cut.height,
    quality.score = best.score,
    diagnostics = diagnostics,
    hclust = data.hclust
  ))
}


#' Compute Binary Partition Diagnostics
#'
#' @keywords internal
compute.binary.partition.diagnostics <- function(data.matrix, 
                                                 partition, 
                                                 hclust,
                                                 cut.height) {
  
  library(cluster)
  
  n <- nrow(data.matrix)
  data.dist <- dist(data.matrix)
  
  # Cluster sizes
  sizes <- table(partition)
  
  # Within-cluster statistics
  within.stats <- lapply(1:2, function(cl) {
    cl.data <- data.matrix[partition == cl, ]
    cl.center <- colMeans(cl.data)
    
    # Within-cluster distances
    if (nrow(cl.data) > 1) {
      cl.dist <- dist(cl.data)
      
      list(
        size = nrow(cl.data),
        mean.dist = mean(cl.dist),
        max.dist = max(cl.dist),
        diameter = max(cl.dist),
        variance = sum(sweep(cl.data, 2, cl.center)^2) / nrow(cl.data)
      )
    } else {
      list(size = 1, mean.dist = 0, max.dist = 0, 
           diameter = 0, variance = 0)
    }
  })
  
  # Between-cluster statistics
  center1 <- colMeans(data.matrix[partition == 1, ])
  center2 <- colMeans(data.matrix[partition == 2, ])
  between.dist <- sqrt(sum((center1 - center2)^2))
  
  # Gap ratio: how much larger is the between-cluster distance
  # compared to within-cluster distances?
  avg.within.diameter <- mean(c(
    within.stats[[1]]$diameter,
    within.stats[[2]]$diameter
  ))
  gap.ratio <- between.dist / avg.within.diameter
  
  # Silhouette
  sil <- silhouette(partition, data.dist)
  avg.sil <- mean(sil[, 3])
  
  # Height analysis
  merge.heights <- hclust$height
  
  # Maximum height within each cluster
  # (find the tallest merge that stays within each cluster)
  max.height.cluster1 <- max(merge.heights[merge.heights < cut.height])
  max.height.cluster2 <- max(merge.heights[merge.heights < cut.height])
  
  # Gap: how much larger is the merge height compared to within-cluster heights?
  height.gap <- cut.height - max(max.height.cluster1, max.height.cluster2)
  height.gap.ratio <- cut.height / max(max.height.cluster1, max.height.cluster2)
  
  # Dunn index: min(between) / max(within)
  # Higher is better
  dunn.index <- min(as.matrix(data.dist)[partition == 1, partition == 2]) /
                max(within.stats[[1]]$diameter, within.stats[[2]]$diameter)
  
  return(list(
    cluster.sizes = sizes,
    within.stats = within.stats,
    between.distance = between.dist,
    gap.ratio = gap.ratio,
    height.gap = height.gap,
    height.gap.ratio = height.gap.ratio,
    avg.silhouette = avg.sil,
    dunn.index = dunn.index,
    silhouette.full = sil
  ))
}


#' Plot Binary Partition Diagnostics
#'
#' @keywords internal
plot.binary.partition.diagnostics <- function(data.matrix,
                                              partition,
                                              hclust,
                                              cut.height,
                                              scores,
                                              test.heights,
                                              diagnostics) {
  
  library(ComplexHeatmap)
  library(circlize)
  
  # Create 2x3 layout
  layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 2, byrow = TRUE))
  
  # 1. Dendrogram with cut
  plot(hclust, labels = FALSE, 
       main = "Dendrogram with Binary Partition",
       xlab = "", sub = "")
  abline(h = cut.height, col = "red", lwd = 3, lty = 1)
  
  # Add height annotations
  max.within.height <- max(hclust$height[hclust$height < cut.height])
  abline(h = max.within.height, col = "blue", lwd = 2, lty = 2)
  
  legend("topright", 
         legend = c("Partition cut", "Max within-cluster height"),
         col = c("red", "blue"), lwd = c(3, 2), lty = c(1, 2),
         bg = "white")
  
  # 2. Height gap analysis
  valid.scores <- scores[!is.na(scores) & scores > 0]
  valid.heights <- test.heights[!is.na(scores) & scores > 0]
  
  if (length(valid.scores) > 0) {
    plot(valid.heights, valid.scores, type = "l", lwd = 2,
         xlab = "Cut Height", ylab = "Partition Quality Score",
         main = "Partition Quality vs Cut Height")
    abline(v = cut.height, col = "red", lwd = 2, lty = 2)
    points(cut.height, max(valid.scores), col = "red", pch = 19, cex = 2)
  }
  
  # 3. Silhouette plot
  plot(diagnostics$silhouette.full, col = c("blue", "red"),
       main = "Silhouette Plot")
  
  # 4. Cluster size comparison
  barplot(diagnostics$cluster.sizes,
          names.arg = c("Cluster 1", "Cluster 2"),
          col = c("blue", "red"),
          main = "Cluster Sizes",
          ylab = "Number of Vertices")
  
  # 5. Diagnostic statistics
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  y.pos <- 0.95
  line.height <- 0.08
  
  text(0.5, y.pos, "Binary Partition Diagnostics", 
       font = 2, cex = 1.2)
  y.pos <- y.pos - line.height * 1.5
  
  stats.text <- c(
    sprintf("Cluster sizes: %d vs %d", 
            diagnostics$cluster.sizes[1],
            diagnostics$cluster.sizes[2]),
    sprintf("Avg silhouette: %.3f", diagnostics$avg.silhouette),
    sprintf("Between-cluster distance: %.3f", diagnostics$between.distance),
    sprintf("Gap ratio (between/within): %.2f", diagnostics$gap.ratio),
    sprintf("Height gap: %.3f", diagnostics$height.gap),
    sprintf("Height gap ratio: %.2f", diagnostics$height.gap.ratio),
    sprintf("Dunn index: %.3f", diagnostics$dunn.index)
  )
  
  for (txt in stats.text) {
    text(0.05, y.pos, txt, pos = 4, cex = 0.9)
    y.pos <- y.pos - line.height
  }
  
  # Interpretation
  y.pos <- y.pos - line.height
  text(0.5, y.pos, "Interpretation:", font = 2, cex = 1.1)
  y.pos <- y.pos - line.height
  
  if (diagnostics$gap.ratio > 3) {
    text(0.05, y.pos, "STRONG binary structure detected!", 
         pos = 4, col = "darkgreen", font = 2)
  } else if (diagnostics$gap.ratio > 1.5) {
    text(0.05, y.pos, "MODERATE binary structure", 
         pos = 4, col = "orange")
  } else {
    text(0.05, y.pos, "WEAK binary structure", 
         pos = 4, col = "red")
  }
  
  layout(1)
  
  # Also create a heatmap with the partition
  row.ha <- rowAnnotation(
    Cluster = factor(partition),
    col = list(Cluster = c("1" = "blue", "2" = "red"))
  )
  
  col.fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  ht <- Heatmap(
    data.matrix[order(partition), ],
    name = "Co-mono",
    col = col.fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    row_split = factor(partition[order(partition)]),
    column_title = "Binary Partition: Data Matrix",
    right_annotation = row.ha,
    show_row_names = FALSE,
    row_title = c("Cluster 1", "Cluster 2"),
    row_title_gp = gpar(col = c("blue", "red"), fontsize = 14)
  )
  
  draw(ht)
}

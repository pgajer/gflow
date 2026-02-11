#' Extract Stable Clusters via Maximal Gap Cluster Peeling
#'
#' @description
#' Implements Maximal Gap Cluster Peeling (MGCP) for identifying stable
#' vertex clusters in co-monotonicity networks.
#'
#' @details
#' In the context of gradient flow analysis, MGCP identifies stable regions
#' of the response surface by iteratively extracting vertex groups with
#' maximal stability (gap size).
#'
#' @param data.matrix Data matrix (vertices x features)
#' @param linkage.method Hierarchical clustering linkage. Default "ward.D2"
#' @param distance.metric Distance metric for dist(). Default "euclidean"
#' @param min.cluster.size Minimum number of elements per extracted cluster
#' @param max.iterations Maximum number of clusters to extract
#' @param min.gap.size Minimum gap size to consider a cluster stable
#' @param min.remaining.size Stop when remaining data has fewer elements
#' @param plot.diagnostics Whether to plot diagnostic figures
#' @param verbose Print progress information
#'
#' @return List containing:
#'   \describe{
#'     \item{clusters}{Integer vector of cluster assignments (0 = unassigned)}
#'     \item{extracted.clusters}{List with details of each extracted cluster}
#'     \item{n.clusters}{Number of clusters extracted}
#'     \item{n.unassigned}{Number of elements not assigned to any cluster}
#'   }
#'
#' @details
#' The algorithm works as follows:
#' 1. Compute hierarchical clustering on current data
#' 2. Identify cluster with largest gap (excluding root node)
#' 3. Extract that cluster if it meets criteria
#' 4. Remove extracted elements from data
#' 5. Repeat on remaining data until stopping criterion
#'
#' This iterative approach finds stable clusters at different heights in the
#' dendrogram, unlike methods that cut at a single level.
#'
#' @export
extract.stable.clusters <- function(data.matrix,
                                    linkage.method = "ward.D2",
                                    distance.metric = "euclidean",
                                    min.cluster.size = 10,
                                    max.iterations = 20,
                                    min.gap.size = NULL,
                                    min.remaining.size = 10,
                                    plot.diagnostics = TRUE,
                                    verbose = TRUE) {
  
  n.total <- nrow(data.matrix)
  
  if (is.null(min.gap.size)) {
    # Auto-detect based on initial clustering
    initial.dist <- dist(data.matrix, method = distance.metric)
    initial.hclust <- hclust(initial.dist, method = linkage.method)
    
    # Use median gap as threshold
    gaps <- compute.cluster.gaps(initial.hclust, min.cluster.size)
    if (nrow(gaps) > 0) {
      min.gap.size <- quantile(gaps$gap.size[gaps$cluster.id != nrow(initial.hclust$merge)], 
                                0.5, na.rm = TRUE)
      if (verbose) {
        cat("Auto-detected min.gap.size:", round(min.gap.size, 3), "\n\n")
      }
    } else {
      min.gap.size <- 0
    }
  }
  
  # Initialize tracking
  cluster.assignments <- integer(n.total)
  names(cluster.assignments) <- rownames(data.matrix)
  if (is.null(names(cluster.assignments))) {
    names(cluster.assignments) <- paste0("V", seq_len(n.total))
  }
  
  extracted.clusters <- list()
  remaining.indices <- seq_len(n.total)
  iteration <- 0
  
  if (verbose) {
    cat("Starting iterative cluster extraction...\n")
    cat("Total elements:", n.total, "\n")
    cat("Min cluster size:", min.cluster.size, "\n")
    cat("Min gap size:", round(min.gap.size, 3), "\n\n")
  }
  
  # Iterative extraction
  while (iteration < max.iterations && length(remaining.indices) >= min.remaining.size) {
    iteration <- iteration + 1
    
    if (verbose) {
      cat("Iteration", iteration, "- Remaining elements:", 
          length(remaining.indices), "\n")
    }
    
    # Cluster remaining data
    current.data <- data.matrix[remaining.indices, , drop = FALSE]
    current.dist <- dist(current.data, method = distance.metric)
    current.hclust <- hclust(current.dist, method = linkage.method)
    
    # Find stable clusters
    gaps <- compute.cluster.gaps(current.hclust, min.cluster.size)
    
    if (nrow(gaps) == 0) {
      if (verbose) {
        cat("  No clusters meeting size requirement\n")
      }
      break
    }
    
    # Exclude root node (last merge)
    n.merges <- nrow(current.hclust$merge)
    gaps <- gaps[gaps$cluster.id != n.merges, ]
    
    if (nrow(gaps) == 0) {
      if (verbose) {
        cat("  Only root cluster available\n")
      }
      break
    }
    
    # Find cluster with largest gap
    best.idx <- which.max(gaps$gap.size)
    best.cluster <- gaps[best.idx, ]
    
    if (verbose) {
      cat("  Best cluster: ID =", best.cluster$cluster.id,
          ", size =", best.cluster$cluster.size,
          ", height =", round(best.cluster$cluster.height, 3),
          ", gap =", round(best.cluster$gap.size, 3), "\n")
    }
    
    # Check if gap meets threshold
    if (best.cluster$gap.size < min.gap.size) {
      if (verbose) {
        cat("  Gap too small (", round(best.cluster$gap.size, 3), 
            " < ", round(min.gap.size, 3), ")\n")
      }
      break
    }
    
    # Extract cluster members
    cluster.members <- extract.cluster.members(current.hclust, 
                                                best.cluster$cluster.id)
    
    # Map back to original indices
    original.indices <- remaining.indices[cluster.members]
    
    # Assign to cluster
    cluster.assignments[original.indices] <- iteration
    
    # Store cluster info
    extracted.clusters[[iteration]] <- list(
      cluster.id = iteration,
      members = original.indices,
      size = length(original.indices),
      height = best.cluster$cluster.height,
      gap.size = best.cluster$gap.size,
      iteration = iteration
    )
    
    # Remove from remaining
    remaining.indices <- remaining.indices[-cluster.members]
    
    if (verbose) {
      cat("  Extracted cluster", iteration, "with", 
          length(original.indices), "elements\n")
      cat("  Remaining:", length(remaining.indices), "elements\n\n")
    }
  }
  
  if (verbose) {
    cat("Extraction complete!\n")
    cat("Total clusters extracted:", length(extracted.clusters), "\n")
    cat("Elements assigned:", sum(cluster.assignments > 0), "\n")
    cat("Elements unassigned:", sum(cluster.assignments == 0), "\n\n")
  }
  
  # Create summary statistics
  summary.stats <- do.call(rbind, lapply(extracted.clusters, function(cl) {
    data.frame(
      cluster.id = cl$cluster.id,
      cluster.size = cl$size,
      cluster.height = cl$height,
      gap.size = cl$gap.size,
      stringsAsFactors = FALSE
    )
  }))
  
  # Diagnostics
  if (plot.diagnostics && length(extracted.clusters) > 0) {
    plot_cluster_extraction_diagnostics(
      data.matrix,
      cluster.assignments,
      extracted.clusters,
      summary.stats
    )
  }
  
  return(list(
    clusters = cluster.assignments,
    extracted.clusters = extracted.clusters,
    summary.stats = summary.stats,
    n.clusters = length(extracted.clusters),
    n.unassigned = sum(cluster.assignments == 0),
    linkage.method = linkage.method,
    distance.metric = distance.metric
  ))
}


#' Compute Cluster Gaps from Hierarchical Clustering
#'
#' Internal function to compute gap sizes for all clusters in dendrogram
#'
#' @keywords internal
compute.cluster.gaps <- function(hclust.obj, min.size = 2) {
  
  merge.matrix <- hclust.obj$merge
  heights <- hclust.obj$height
  n.merges <- nrow(merge.matrix)
  n.elements <- n.merges + 1
  
  # Initialize tracking
  cluster.list <- vector("list", n.merges)
  cluster.birth.height <- numeric(n.merges)
  cluster.death.height <- numeric(n.merges)
  
  # Process each merge
  for (i in seq_len(n.merges)) {
    left <- merge.matrix[i, 1]
    right <- merge.matrix[i, 2]
    current.height <- heights[i]
    
    # Get left members
    if (left < 0) {
      left.members <- abs(left)
    } else {
      left.members <- cluster.list[[left]]
      cluster.death.height[left] <- current.height
    }
    
    # Get right members
    if (right < 0) {
      right.members <- abs(right)
    } else {
      right.members <- cluster.list[[right]]
      cluster.death.height[right] <- current.height
    }
    
    # Combine
    cluster.list[[i]] <- c(left.members, right.members)
    cluster.birth.height[i] <- current.height
  }
  
  # Root never dies
  cluster.death.height[n.merges] <- max(heights) * 2
  
  # Compute gaps
  gap.sizes <- cluster.death.height - cluster.birth.height
  cluster.sizes <- sapply(cluster.list, length)
  
  # Create stats data frame
  stats <- data.frame(
    cluster.id = seq_len(n.merges),
    cluster.size = cluster.sizes,
    cluster.height = cluster.birth.height,
    gap.size = gap.sizes,
    stringsAsFactors = FALSE
  )
  
  # Filter by size
  stats <- stats[stats$cluster.size >= min.size, ]
  
  # Sort by gap (descending)
  stats <- stats[order(-stats$gap.size), ]
  rownames(stats) <- NULL
  
  return(stats)
}


#' Extract Cluster Members from Hierarchical Clustering
#'
#' Internal function to get all leaf members of a cluster
#'
#' @keywords internal
extract.cluster.members <- function(hclust.obj, cluster.id) {
  
  merge.matrix <- hclust.obj$merge
  n.merges <- nrow(merge.matrix)
  
  # Build cluster membership recursively
  get.members <- function(node.id) {
    if (node.id < 0) {
      # Leaf node
      return(abs(node.id))
    } else {
      # Internal node - recurse
      left <- merge.matrix[node.id, 1]
      right <- merge.matrix[node.id, 2]
      return(c(get.members(left), get.members(right)))
    }
  }
  
  members <- get.members(cluster.id)
  return(sort(members))
}


#' Plot Cluster Extraction Diagnostics
#'
#' Creates diagnostic visualizations for iterative cluster extraction
#'
#' @keywords internal
plot_cluster_extraction_diagnostics <- function(data.matrix,
                                                cluster.assignments,
                                                extracted.clusters,
                                                summary.stats) {
  n.clusters <- length(extracted.clusters)
  
  # Create 2x2 layout
  par(mfrow = c(2, 2))
  
  # 1. Gap sizes by extraction order
  plot(summary.stats$cluster.id, summary.stats$gap.size,
       type = "b", pch = 19, col = "blue", lwd = 2,
       xlab = "Extraction Order", ylab = "Gap Size",
       main = "Cluster Stability by Extraction Order")
  grid()
  
  # 2. Cluster sizes
  barplot(summary.stats$cluster.size,
          names.arg = summary.stats$cluster.id,
          col = rainbow(n.clusters),
          xlab = "Cluster ID", ylab = "Size",
          main = "Cluster Sizes")
  
  # 3. Height vs Gap
  plot(summary.stats$cluster.height, summary.stats$gap.size,
       pch = 19, cex = 1.5, col = rainbow(n.clusters),
       xlab = "Cluster Height", ylab = "Gap Size",
       main = "Cluster Height vs Stability")
  text(summary.stats$cluster.height, summary.stats$gap.size,
       labels = summary.stats$cluster.id, pos = 3, cex = 0.8)
  grid()
  
  # 4. Summary text
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  y.pos <- 0.95
  line.height <- 0.08
  
  text(0.5, y.pos, "Extraction Summary", font = 2, cex = 1.3)
  y.pos <- y.pos - line.height * 1.5
  
  summary.text <- c(
    sprintf("Total elements: %d", length(cluster.assignments)),
    sprintf("Clusters extracted: %d", n.clusters),
    sprintf("Elements assigned: %d", sum(cluster.assignments > 0)),
    sprintf("Elements unassigned: %d", sum(cluster.assignments == 0)),
    "",
    sprintf("Largest cluster: %d elements", max(summary.stats$cluster.size)),
    sprintf("Smallest cluster: %d elements", min(summary.stats$cluster.size)),
    sprintf("Max gap: %.3f", max(summary.stats$gap.size)),
    sprintf("Min gap: %.3f", min(summary.stats$gap.size))
  )
  
  for (txt in summary.text) {
    if (txt != "") {
      text(0.05, y.pos, txt, pos = 4, cex = 0.9)
    }
    y.pos <- y.pos - line.height
  }
  
  par(mfrow = c(1, 1))
  
  # Create heatmap if clusters found
  if (n.clusters > 0) {
    # Order by cluster assignment
    assigned <- cluster.assignments > 0
    if (sum(assigned) > 0) {
      order.idx <- order(cluster.assignments)
      
      row.ha <- ComplexHeatmap::rowAnnotation(
        Cluster = factor(cluster.assignments[order.idx],
                        levels = 0:n.clusters),
        col = list(Cluster = setNames(
          c("gray", rainbow(n.clusters)),
          as.character(0:n.clusters)
        )),
        show_legend = TRUE,
        annotation_name_side = "top"
      )
      
      col.fun <- circlize::colorRamp2(
        seq(min(data.matrix), max(data.matrix), length.out = 5),
        c("blue", "lightblue", "white", "pink", "red")
      )
      
      ht <- ComplexHeatmap::Heatmap(
        data.matrix[order.idx, ],
        name = "Value",
        col = col.fun,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        row_split = factor(cluster.assignments[order.idx],
                          levels = 0:n.clusters),
        column_title = "Extracted Clusters",
        right_annotation = row.ha,
        show_row_names = FALSE,
        row_title = c("Unassigned", paste("Cluster", 1:n.clusters))
      )
      
      ComplexHeatmap::draw(ht)
    }
  }
}


#' Analyze Extracted Cluster Characteristics
#'
#' Computes detailed statistics for each extracted cluster
#'
#' @param data.matrix Original data matrix
#' @param extraction.result Result from extract.stable.clusters()
#'
#' @return Data frame with cluster characteristics
#'
#' @export
analyze.extracted.clusters <- function(data.matrix, extraction.result) {
  clusters <- extraction.result$clusters
  n.clusters <- extraction.result$n.clusters
  
  if (n.clusters == 0) {
    warning("No clusters to analyze")
    return(NULL)
  }
  
  # Compute distance matrix once
  data.dist <- dist(data.matrix)
  
  # Analyze each cluster
  results <- list()
  
  for (cl in seq_len(n.clusters)) {
    cl.members <- which(clusters == cl)
    cl.data <- data.matrix[cl.members, , drop = FALSE]
    
    # Within-cluster statistics
    if (nrow(cl.data) > 1) {
      cl.dist <- dist(cl.data)
      cl.center <- colMeans(cl.data)
      
      # Coherence metrics
      cl.cors <- cor(t(cl.data))
      diag(cl.cors) <- NA
      
      # Safely compute median correlation
      cors.vec <- as.vector(cl.cors[upper.tri(cl.cors)])
      cors.vec <- cors.vec[!is.na(cors.vec)]
      
      if (length(cors.vec) > 0) {
        median.cor <- median(abs(cors.vec))
      } else {
        median.cor <- NA_real_
      }
      
      # SVD for rank structure
      svd.cl <- svd(cl.data)
      var.exp <- svd.cl$d^2 / sum(svd.cl$d^2)
      
      results[[cl]] <- list(
        cluster.id = cl,
        size = nrow(cl.data),
        mean.within.dist = mean(cl.dist),
        median.correlation = median.cor,
        rank1.variance = var.exp[1],
        rank2.variance = sum(var.exp[1:min(2, length(var.exp))]),
        effective.rank = sum(var.exp > 0.01)
      )
    } else {
      results[[cl]] <- list(
        cluster.id = cl,
        size = 1,
        mean.within.dist = 0,
        median.correlation = NA_real_,
        rank1.variance = 1,
        rank2.variance = 1,
        effective.rank = 1
      )
    }
  }
  
  # Convert to data frame
  if (length(results) == 0) {
    return(data.frame(
      cluster.id = integer(0),
      size = integer(0),
      mean.within.dist = numeric(0),
      median.correlation = numeric(0),
      rank1.variance = numeric(0),
      rank2.variance = numeric(0),
      effective.rank = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract each component separately to ensure proper lengths
  cluster.ids <- sapply(results, function(x) x$cluster.id)
  sizes <- sapply(results, function(x) x$size)
  mean.dists <- sapply(results, function(x) x$mean.within.dist)
  median.cors <- sapply(results, function(x) x$median.correlation)
  rank1.vars <- sapply(results, function(x) x$rank1.variance)
  rank2.vars <- sapply(results, function(x) x$rank2.variance)
  eff.ranks <- sapply(results, function(x) x$effective.rank)
  
  results.df <- data.frame(
    cluster.id = cluster.ids,
    size = sizes,
    mean.within.dist = mean.dists,
    median.correlation = median.cors,
    rank1.variance = rank1.vars,
    rank2.variance = rank2.vars,
    effective.rank = eff.ranks,
    stringsAsFactors = FALSE
  )
  
  return(results.df)
}

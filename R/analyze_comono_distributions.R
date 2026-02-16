#' Analyze Co-Monotonicity Coefficient Distributions by Cell
#'
#' Examines the distribution of co-monotonicity coefficients within each super-cell
#' to identify patterns in the strength and direction of associations.
#'
#' @param cm Co-monotonicity matrix (samples x phylotypes).
#' @param partition Vector of super-cell membership.
#' @param cell.ids Optional vector of specific cells to analyze. If NULL, analyzes all.
#' @param response.values Optional vector of response values (e.g., sPTB prevalence) for
#'   each sample, used to understand heterogeneity within cells.
#'
#' @return A list containing diagnostic information for each cell.
#'
analyze.comono.distributions <- function(cm,
                                         partition,
                                         cell.ids = NULL,
                                         response.values = NULL) {

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

  if (is.data.frame(cm)) {
    cm <- as.matrix(cm)
  }

  ## Get cell IDs
  if (is.null(cell.ids)) {
    cell.ids <- sort(unique(partition))
  }

  ## Initialize results list
  results <- list()

  for (cell.id in cell.ids) {
    ## Extract cell data
    cell.rows <- which(partition == cell.id)
    cell.matrix <- cm[cell.rows, , drop = FALSE]
    n.cell.samples <- nrow(cell.matrix)
    n.cell.phylotypes <- ncol(cell.matrix)

    ## Create sample IDs if row names don't exist
    if (is.null(rownames(cell.matrix))) {
      sample.ids <- paste0("Sample", seq_len(n.cell.samples))
    } else {
      sample.ids <- rownames(cell.matrix)
    }

    ## Create phylotype IDs if column names don't exist
    if (is.null(colnames(cell.matrix))) {
      phylotype.ids <- paste0("Phylotype", seq_len(n.cell.phylotypes))
    } else {
      phylotype.ids <- colnames(cell.matrix)
    }

    ## Compute phylotype-level statistics
    phylotype.stats <- data.frame(
      phylotype = phylotype.ids,
      mean.comono = colMeans(cell.matrix, na.rm = TRUE),
      median.comono = apply(cell.matrix, 2, median, na.rm = TRUE),
      sd.comono = apply(cell.matrix, 2, sd, na.rm = TRUE),
      min.comono = apply(cell.matrix, 2, min, na.rm = TRUE),
      max.comono = apply(cell.matrix, 2, max, na.rm = TRUE),
      range.comono = apply(cell.matrix, 2, function(x) diff(range(x, na.rm = TRUE))),
      prop.near.zero = apply(cell.matrix, 2, function(x) mean(abs(x) < 0.1, na.rm = TRUE)),
      prop.extreme = apply(cell.matrix, 2, function(x) mean(abs(x) > 0.8, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    rownames(phylotype.stats) <- NULL

    ## Compute sample-level statistics
    sample.stats <- data.frame(
      sample.id = sample.ids,
      mean.abs.comono = rowMeans(abs(cell.matrix), na.rm = TRUE),
      median.abs.comono = apply(abs(cell.matrix), 1, median, na.rm = TRUE),
      prop.near.zero = apply(cell.matrix, 1, function(x) mean(abs(x) < 0.1, na.rm = TRUE)),
      prop.extreme = apply(cell.matrix, 1, function(x) mean(abs(x) > 0.8, na.rm = TRUE)),
      prop.positive = apply(cell.matrix, 1, function(x) mean(x > 0.1, na.rm = TRUE)),
      prop.negative = apply(cell.matrix, 1, function(x) mean(x < -0.1, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    rownames(sample.stats) <- NULL

    ## Add response values if provided
    if (!is.null(response.values)) {
      sample.stats$response <- response.values[cell.rows]
    }

    ## Overall matrix statistics
    matrix.stats <- list(
      n.samples = n.cell.samples,
      n.phylotypes = n.cell.phylotypes,
      overall.mean = mean(cell.matrix, na.rm = TRUE),
      overall.median = median(cell.matrix, na.rm = TRUE),
      overall.sd = sd(as.vector(cell.matrix), na.rm = TRUE),
      prop.near.zero = mean(abs(cell.matrix) < 0.1, na.rm = TRUE),
      prop.extreme = mean(abs(cell.matrix) > 0.8, na.rm = TRUE),
      prop.positive = mean(cell.matrix > 0.1, na.rm = TRUE),
      prop.negative = mean(cell.matrix < -0.1, na.rm = TRUE)
    )

    ## Store results
    results[[as.character(cell.id)]] <- list(
      phylotype.stats = phylotype.stats,
      sample.stats = sample.stats,
      matrix.stats = matrix.stats,
      cell.matrix = cell.matrix
    )
  }

  class(results) <- c("comono_analysis", "list")
  return(results)
}

#' Create Diagnostic Plots for Co-Monotonicity Distribution Analysis
#'
#' @param x Output from analyze.comono.distributions().
#' @param cell.id Specific cell ID to plot.
#' @param ... Additional arguments (currently ignored).
#'
#' @export
plot.comono.diagnostics <- function(x, cell.id, ...) {
  analysis.results <- x

  if (!inherits(analysis.results, "comono_analysis")) {
    stop("analysis.results must be output from analyze.comono.distributions()")
  }

  cell.id.char <- as.character(cell.id)
  if (!cell.id.char %in% names(analysis.results)) {
    stop("Cell ID not found in analysis results")
  }

  cell.data <- analysis.results[[cell.id.char]]

  ## Set up multi-panel plot
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), mgp = c(2.5, 0.7, 0))

  ## 1. Overall distribution of co-monotonicity values
  hist(cell.data$cell.matrix,
       breaks = 50,
       col = "lightblue",
       border = "white",
       main = paste("Cell", cell.id, "- Overall Distribution"),
       xlab = "Co-Monotonicity Coefficient",
       ylab = "Frequency")
  abline(v = c(-0.1, 0.1), col = "red", lty = 2, lwd = 2)

  ## 2. Distribution of phylotype means
  hist(cell.data$phylotype.stats$mean.comono,
       breaks = 30,
       col = "lightgreen",
       border = "white",
       main = "Phylotype Mean Co-Monotonicity",
       xlab = "Mean Co-Monotonicity",
       ylab = "Number of Phylotypes")
  abline(v = 0, col = "red", lty = 2, lwd = 2)

  ## 3. Proportion of near-zero values per phylotype
  hist(cell.data$phylotype.stats$prop.near.zero,
       breaks = 30,
       col = "lightyellow",
       border = "white",
       main = "Near-Zero Proportions (Phylotypes)",
       xlab = "Proportion |comono| < 0.1",
       ylab = "Number of Phylotypes")

  ## 4. Proportion of extreme values per phylotype
  hist(cell.data$phylotype.stats$prop.extreme,
       breaks = 30,
       col = "lightcoral",
       border = "white",
       main = "Extreme Value Proportions (Phylotypes)",
       xlab = "Proportion |comono| > 0.8",
       ylab = "Number of Phylotypes")

  ## 5. Sample heterogeneity: proportion near zero per sample
  plot(cell.data$sample.stats$prop.near.zero,
       cell.data$sample.stats$prop.extreme,
       pch = 19,
       col = rgb(0, 0, 1, 0.5),
       main = "Sample Heterogeneity",
       xlab = "Proportion Near Zero",
       ylab = "Proportion Extreme")
  abline(h = 0.5, v = 0.5, col = "red", lty = 2)

  ## 6. Response vs proportion extreme (if response available)
  if ("response" %in% names(cell.data$sample.stats)) {
    plot(cell.data$sample.stats$response,
         cell.data$sample.stats$prop.extreme,
         pch = 19,
         col = rgb(0, 0, 1, 0.5),
         main = "Response vs Extremity",
         xlab = "Response Value",
         ylab = "Proportion Extreme")

    ## Add smoothing line
    lo <- loess(prop.extreme ~ response, data = cell.data$sample.stats)
    response.order <- order(cell.data$sample.stats$response)
    lines(cell.data$sample.stats$response[response.order],
          predict(lo)[response.order],
          col = "red",
          lwd = 2)
  } else {
    ## Alternative: show relationship between positive and negative
    plot(cell.data$sample.stats$prop.positive,
         cell.data$sample.stats$prop.negative,
         pch = 19,
         col = rgb(0, 0, 1, 0.5),
         main = "Directional Balance",
         xlab = "Proportion Positive",
         ylab = "Proportion Negative")
    abline(a = 1, b = -1, col = "red", lty = 2)
  }

  par(mfrow = c(1, 1))
}

#' Identify and Save High-Discrepancy Vertex-Phylotype Pairs for Debugging
#'
#' Finds pairs where standard and sign co-monotonicity coefficients differ most
#' and saves them to a CSV file for C++ debugging.
#'
#' @param cm.standard Standard co-monotonicity matrix (vertices x phylotypes)
#' @param cm.sign Sign co-monotonicity matrix (vertices x phylotypes)
#' @param partition Partition vector
#' @param cell.id Cell to analyze
#' @param n.pairs Number of pairs to return
#' @param debug.file Path to save debug pairs CSV. Default uses fixed debug directory.
#'
#' @return Data frame with high-discrepancy pairs
find.discrepant.pairs <- function(cm.standard,
                                   cm.sign,
                                   partition = NULL,
                                   cell.id = NULL,
                                   n.pairs = 10,
                                   debug.file = NULL) {

  ## Use fixed debug directory if not specified
  if (is.null(debug.file)) {
    debug.dir <- "/tmp/comono_debugging_dir"

    ## Create directory if it doesn't exist
    if (!dir.exists(debug.dir)) {
      dir.create(debug.dir, recursive = TRUE)
      cat("Created debug directory:", debug.dir, "\n")
    }

    debug.file <- file.path(debug.dir, "comono_debug_pairs.csv")
  }

  ## Determine which rows to analyze
  if (!is.null(partition) && !is.null(cell.id)) {
    cell.rows <- which(partition == cell.id)
    cat("Analyzing cell", cell.id, "with", length(cell.rows), "samples\n")
  } else {
    cell.rows <- seq_len(nrow(cm.standard))
    cat("Analyzing all", length(cell.rows), "samples\n")
  }

  ## Compute absolute difference matrix for selected rows
  diff.matrix <- abs(cm.standard[cell.rows, ] - cm.sign[cell.rows, ])

  ## Find pairs with largest discrepancies
  discrepancies <- list()

  for (i in seq_along(cell.rows)) {
    for (j in seq_len(ncol(cm.standard))) {
      discrepancies[[length(discrepancies) + 1]] <- list(
        vertex.idx.local = i,
        phylotype.idx = j,
        vertex.idx.global = cell.rows[i],
        phylotype.name = if (!is.null(colnames(cm.standard))) {
          colnames(cm.standard)[j]
        } else {
          paste0("Phylotype", j)
        },
        standard.value = cm.standard[cell.rows[i], j],
        sign.value = cm.sign[cell.rows[i], j],
        difference = diff.matrix[i, j]
      )
    }
  }

  ## Convert to data frame
  discrepancies.df <- do.call(rbind, lapply(discrepancies, as.data.frame))

  ## Sort by difference and return top pairs
  discrepancies.df <- discrepancies.df[order(-discrepancies.df$difference), ]
  top.pairs <- head(discrepancies.df, n.pairs)

  ## Add case categories for interpretation
  top.pairs$case <- with(top.pairs, {
    ifelse(
      abs(standard.value) > 0.8 & abs(sign.value) < 0.3,
      "extreme_standard_moderate_sign",
      ifelse(
        abs(standard.value) < 0.3 & abs(sign.value) > 0.8,
        "moderate_standard_extreme_sign",
        ifelse(
          sign(standard.value) != sign(sign.value),
          "opposite_signs",
          "both_extreme_same_direction"
        )
      )
    )
  })

  ## Save to CSV for C++ consumption
  ## Use 0-based indexing for C++
  debug.csv <- data.frame(
    vertex_idx = top.pairs$vertex.idx.global - 1,  # Convert to 0-based
    phylotype_idx = top.pairs$phylotype.idx - 1,   # Convert to 0-based
    phylotype_name = top.pairs$phylotype.name,
    standard_value = top.pairs$standard.value,
    sign_value = top.pairs$sign.value,
    difference = top.pairs$difference,
    case = top.pairs$case
  )

  write.csv(debug.csv, debug.file, row.names = FALSE)
  cat("Debug pairs saved to:", debug.file, "\n")

  ## Also save the R-friendly version with 1-based indexing
  debug.file.r <- sub("\\.csv$", "_R.csv", debug.file)
  write.csv(top.pairs, debug.file.r, row.names = FALSE)
  cat("R version (1-based indices) saved to:", debug.file.r, "\n")

  ## Print summary
  cat("\nTop", min(5, nrow(top.pairs)), "discrepancies:\n")
  print(top.pairs[1:min(5, nrow(top.pairs)),
                   c("vertex.idx.global", "phylotype.name",
                     "standard.value", "sign.value", "difference", "case")])

  return(top.pairs)
}

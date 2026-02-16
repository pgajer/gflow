#' Compute Partition Cell Response Variable Statistics
#'
#' Computes summary statistics for response values within each cell of a
#' partition.
#'
#' @param partition Numeric vector of cluster membership assignments (1-based indexing).
#'   Each element indicates which partition cell the corresponding vertex belongs to.
#' @param y Numeric vector of response values, same length as partition.
#' @param compute.mad Logical indicating whether to compute MAD (median absolute
#'   deviation). Default is TRUE.
#'
#' @return A data.frame with one row per super-cell containing:
#'   \item{cell.id}{Super-cell identifier}
#'   \item{n.vertices}{Number of vertices in the super-cell}
#'   \item{min.value}{Minimum response value}
#'   \item{max.value}{Maximum response value}
#'   \item{mean.value}{Mean response value}
#'   \item{median.value}{Median response value}
#'   \item{range}{Range (max - min)}
#'   \item{iqr}{Interquartile range}
#'   \item{mad}{Median absolute deviation (if compute.mad = TRUE)}
#'
#' @details
#' This function summarizes the distribution of response values within each
#' super-cell of a partition. The statistics provide insight into the variability
#' and central tendency of the response across the geometric decomposition.
#'
compute.partition.cell.stats <- function(partition, y, compute.mad = TRUE) {

  ## Input validation
  if (!is.numeric(partition)) {
    stop("partition must be a numeric vector")
  }
  if (!is.numeric(y)) {
    stop("y must be a numeric vector")
  }
  if (length(partition) != length(y)) {
    stop("partition and y must have the same length")
  }
  if (any(is.na(partition))) {
    stop("partition cannot contain NA values")
  }

  ## Get unique super-cell IDs
  cell.ids <- sort(unique(partition))
  n.cells <- length(cell.ids)

  ## Initialize result data frame
  result <- data.frame(
    cell.id = cell.ids,
    n.vertices = integer(n.cells),
    min.value = numeric(n.cells),
    max.value = numeric(n.cells),
    mean.value = numeric(n.cells),
    median.value = numeric(n.cells),
    range = numeric(n.cells),
    iqr = numeric(n.cells),
    stringsAsFactors = FALSE
  )

  if (compute.mad) {
    result$mad <- numeric(n.cells)
  }

  ## Compute statistics for each super-cell
  for (i in seq_along(cell.ids)) {
    cell.id <- cell.ids[i]
    cell.mask <- partition == cell.id
    cell.values <- y[cell.mask]

    ## Remove NA values for statistical computation
    cell.values.clean <- cell.values[!is.na(cell.values)]

    if (length(cell.values.clean) == 0) {
      warning(sprintf("Super-cell %d has no non-NA values", cell.id))
      result$n.vertices[i] <- sum(cell.mask)
      result$min.value[i] <- NA
      result$max.value[i] <- NA
      result$mean.value[i] <- NA
      result$median.value[i] <- NA
      result$range[i] <- NA
      result$iqr[i] <- NA
      if (compute.mad) {
        result$mad[i] <- NA
      }
      next
    }

    ## Compute summary statistics
    result$n.vertices[i] <- sum(cell.mask)
    result$min.value[i] <- min(cell.values.clean)
    result$max.value[i] <- max(cell.values.clean)
    result$mean.value[i] <- mean(cell.values.clean)
    result$median.value[i] <- median(cell.values.clean)
    result$range[i] <- result$max.value[i] - result$min.value[i]

    ## IQR (handling cases with too few observations)
    if (length(cell.values.clean) >= 4) {
      result$iqr[i] <- IQR(cell.values.clean, na.rm = TRUE)
    } else {
      result$iqr[i] <- NA
    }

    ## MAD if requested
    if (compute.mad) {
      if (length(cell.values.clean) >= 2) {
        result$mad[i] <- mad(cell.values.clean, na.rm = TRUE)
      } else {
        result$mad[i] <- NA
      }
    }
  }

  ## Add informative row names
  rownames(result) <- paste0("cell.", result$cell.id)

  return(result)
}

#' Extract Response Values by Super-Cell
#'
#' Extracts the response values for each super-cell in a partition, returning
#' a list where each element contains all response values belonging to that cell.
#'
#' @param partition Numeric vector of cluster membership assignments (1-based indexing).
#'   Each element indicates which super-cell the corresponding vertex belongs to.
#' @param y Numeric vector of response values, same length as partition.
#'
#' @return A named list with one element per super-cell. Each element is a numeric
#'   vector containing all response values for vertices in that super-cell. List
#'   names correspond to the super-cell identifiers.
#'
#' @details
#' This function provides access to the raw response values within each super-cell
#' of a partition, enabling custom statistical analyses beyond summary statistics.
#' Unlike \code{compute.super.cell.stats}, which returns aggregated summaries,
#' this function returns the complete vector of values for each cell. The result
#' is useful for visualization, distribution analysis, or computing custom statistics
#' that require access to individual observations.
#'
#' The function preserves NA values in the response vector, allowing the caller
#' to decide how to handle missing data. Super-cells are identified by their
#' unique values in the partition vector, and the returned list is indexed by
#' these identifiers.
#'
#' @examples
#' \dontrun{
#' ## Extract values for each super-cell
#' cell.values <- cell.y.values(partition = sptb.phylo.super.cells,
#'                               y = rel.sptb.hat.winsorized)
#'
#' ## Examine distribution in a specific cell
#' hist(cell.values[[3]], main = "Cell 3 Response Distribution")
#'
#' ## Compute custom statistic (e.g., coefficient of variation)
#' cv.by.cell <- sapply(cell.values, function(vals) {
#'   sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE)
#' })
#'
#' ## Compare distributions across cells
#' boxplot(cell.values, xlab = "Super-Cell", ylab = "Response Value")
#' }
#'
#' @seealso \code{compute.super.cell.stats} for computing summary statistics
#'   rather than extracting raw values.
#'
#' @export
cell.y.values <- function(partition, y) {
    ## Input validation
    if (!is.numeric(partition)) {
        stop("partition must be a numeric vector")
    }
    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }
    if (length(partition) != length(y)) {
        stop("partition and y must have the same length")
    }
    if (any(is.na(partition))) {
        stop("partition cannot contain NA values")
    }
    ## Get unique super-cell IDs
    cell.ids <- sort(unique(partition))
    n.cells <- length(cell.ids)
    ## Initialize result list
    result <- list()
    ## Extract values for each super-cell
    for (i in seq_along(cell.ids)) {
        cell.id <- cell.ids[i]
        cell.mask <- partition == cell.id
        cell.values <- y[cell.mask]
        result[[cell.id]] <- cell.values
    }
    return(result)
}

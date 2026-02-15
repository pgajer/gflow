#' Extract simplex vertex lists from riem_dcx object
#'
#' @param dcx A riem_dcx object
#' @param dim Integer dimension (0 for vertices, 1 for edges, etc.)
#' @return List of integer vectors, each representing a simplex
get_simplices <- function(dcx, dim) {
  # This assumes you have an R interface to extract simplices
  # You may need to implement this in your C++ code
  .Call("S_get_simplices", dcx, as.integer(dim), PACKAGE = "gflow")
}

#' Get metric diagonal values for a dimension
#'
#' @param dcx A riem_dcx object
#' @param dim Integer dimension
#' @return Numeric vector of diagonal metric values
get_metric_diagonal <- function(dcx, dim) {
  .Call("S_get_metric_diagonal", dcx, as.integer(dim), PACKAGE = "gflow")
}

#' Check if simplices are sorted
#'
#' @param simplices List of integer vectors
#' @return Logical TRUE if all are sorted
simplices_are_sorted <- function(simplices) {
  all(vapply(simplices, function(s) {
    all(s[-length(s)] < s[-1])
  }, logical(1)))
}

#' Check for duplicate simplices
#'
#' @param simplices List of integer vectors
#' @return List with has_duplicates (logical) and duplicates (list of duplicated simplices)
check_duplicates <- function(simplices) {
  # Convert to strings for comparison
  simplex_strings <- vapply(simplices, function(s) {
    paste(sort(s), collapse = ",")
  }, character(1))

  duplicated_mask <- duplicated(simplex_strings)
  has_duplicates <- any(duplicated_mask)

  if (has_duplicates) {
    duplicate_indices <- which(simplex_strings %in% simplex_strings[duplicated_mask])
    duplicate_simplices <- simplices[duplicate_indices]
  } else {
    duplicate_simplices <- list()
  }

  list(
    has_duplicates = has_duplicates,
    duplicates = duplicate_simplices,
    duplicate_indices = if (has_duplicates) duplicate_indices else integer(0)
  )
}

#' Verify all faces of a simplex exist in lower dimension
#'
#' @param simplex Integer vector of vertex indices
#' @param lower_simplices List of simplices in dimension p-1
#' @return Logical TRUE if all faces exist
all_faces_exist <- function(simplex, lower_simplices) {
  p <- length(simplex)
  if (p <= 1) return(TRUE)

  # Generate all (p-1)-faces by omitting each vertex
  faces <- lapply(seq_len(p), function(i) {
    sort(simplex[-i])
  })

  # Convert lower_simplices to comparable form
  lower_sorted <- lapply(lower_simplices, sort)

  # Check each face exists
  all(vapply(faces, function(face) {
    any(vapply(lower_sorted, function(ls) {
      identical(face, ls)
    }, logical(1)))
  }, logical(1)))
}

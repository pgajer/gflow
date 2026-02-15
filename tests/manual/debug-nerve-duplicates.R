# Manual debugging script for nerve construction
# Not run during R CMD check or devtools::test()

library(gflow)

# Enable detailed output
options(warn = 1)

# ==============================================================================
# Test 1: Small example with known structure
# ==============================================================================

cat("\n=== Test 1: Small example ===\n")
set.seed(12345)
n <- 30
X <- cbind(rnorm(n), rnorm(n))
y <- rnorm(n)

dcx <- build.nerve.from.knn(
  X = X,
  y = y,
  k = 10,
  max.p = 3,
  use.counting.measure = TRUE,
  directed.knn = FALSE
)

summary_info <- summary(dcx)
print(summary_info)

# Check for duplicates manually
check_duplicates_manual <- function(dcx, dim) {
  simplices <- .Call("S_get_simplices", dcx, as.integer(dim), PACKAGE = "gflow")

  cat(sprintf("\nDimension %d: %d simplices\n", dim, length(simplices)))

  if (length(simplices) == 0) {
    cat("  (empty)\n")
    return(invisible(NULL))
  }

  # Convert to strings for comparison
  simplex_strings <- vapply(simplices, function(s) {
    paste(sort(s), collapse = ",")
  }, character(1))

  # Find duplicates
  dup_mask <- duplicated(simplex_strings)

  if (any(dup_mask)) {
    cat(sprintf("  WARNING: Found %d duplicate simplices!\n", sum(dup_mask)))

    # Show details
    dup_indices <- which(simplex_strings %in% simplex_strings[dup_mask])
    for (idx in dup_indices) {
      cat(sprintf("    Simplex %d: {%s}\n",
                  idx,
                  paste(simplices[[idx]], collapse = ", ")))
    }

    # Show unique duplicated patterns
    unique_dups <- unique(simplex_strings[dup_mask])
    cat(sprintf("\n  Unique duplicated patterns (%d):\n", length(unique_dups)))
    for (pattern in unique_dups) {
      indices <- which(simplex_strings == pattern)
      cat(sprintf("    {%s} appears at indices: %s\n",
                  pattern,
                  paste(indices, collapse = ", ")))
    }
  } else {
    cat("  ✓ No duplicates found\n")
  }

  invisible(list(
    has_duplicates = any(dup_mask),
    n_duplicates = sum(dup_mask),
    duplicate_indices = if (any(dup_mask)) which(dup_mask) else integer(0)
  ))
}

# Check all dimensions
cat("\n=== Checking for duplicates ===\n")
for (p in 0:summary_info$max_dimension) {
  check_duplicates_manual(dcx, p)
}

# ==============================================================================
# Test 2: Verify metric values
# ==============================================================================

cat("\n\n=== Test 2: Metric values ===\n")
for (p in 0:summary_info$max_dimension) {
  metric_diag <- .Call("S_get_metric_diagonal", dcx, as.integer(p),
                       PACKAGE = "gflow")

  cat(sprintf("\nDimension %d:\n", p))
  cat(sprintf("  Min metric: %.6f\n", min(metric_diag)))
  cat(sprintf("  Max metric: %.6f\n", max(metric_diag)))
  cat(sprintf("  Mean metric: %.6f\n", mean(metric_diag)))
  cat(sprintf("  Any negative? %s\n", any(metric_diag < 0)))
  cat(sprintf("  Any infinite? %s\n", any(!is.finite(metric_diag))))

  if (length(metric_diag) <= 10) {
    cat("  All values:", paste(sprintf("%.4f", metric_diag), collapse = ", "), "\n")
  }
}

# ==============================================================================
# Test 3: Face existence check
# ==============================================================================

cat("\n\n=== Test 3: Face existence ===\n")

check_faces <- function(dcx, dim) {
  if (dim == 0) {
    cat("Dimension 0: vertices have no faces\n")
    return(invisible(TRUE))
  }

  simplices_p <- .Call("S_get_simplices", dcx, as.integer(dim),
                       PACKAGE = "gflow")
  simplices_p_minus_1 <- .Call("S_get_simplices", dcx, as.integer(dim - 1),
                                PACKAGE = "gflow")

  cat(sprintf("\nDimension %d: checking %d simplices\n",
              dim, length(simplices_p)))

  if (length(simplices_p) == 0) {
    cat("  (no simplices)\n")
    return(invisible(TRUE))
  }

  all_ok <- TRUE
  missing_faces_count <- 0

  for (i in seq_along(simplices_p)) {
    simplex <- simplices_p[[i]]
    p <- length(simplex)

    # Generate all (p-1)-faces
    faces <- lapply(seq_len(p), function(j) {
      sort(simplex[-j])
    })

    # Check each face exists
    for (face in faces) {
      face_exists <- any(vapply(simplices_p_minus_1, function(s) {
        identical(sort(s), face)
      }, logical(1)))

      if (!face_exists) {
        if (all_ok) {
          cat("  ERROR: Missing faces detected!\n")
          all_ok <- FALSE
        }
        missing_faces_count <- missing_faces_count + 1
        if (missing_faces_count <= 5) {  # Show first 5
          cat(sprintf("    Simplex %d {%s} missing face {%s}\n",
                      i,
                      paste(simplex, collapse = ", "),
                      paste(face, collapse = ", ")))
        }
      }
    }
  }

  if (all_ok) {
    cat("  ✓ All faces exist\n")
  } else {
    cat(sprintf("  Total missing faces: %d\n", missing_faces_count))
  }

  invisible(all_ok)
}

for (p in 1:summary_info$max_dimension) {
  check_faces(dcx, p)
}

# ==============================================================================
# Test 4: Clustered data (more likely to have high-dimensional simplices)
# ==============================================================================

cat("\n\n=== Test 4: Clustered data ===\n")
set.seed(42)
n <- 100
# Three tight clusters
centers <- rbind(c(0, 0), c(5, 5), c(0, 5))
cluster_id <- sample(1:3, n, replace = TRUE, prob = c(0.4, 0.3, 0.3))
X_cluster <- t(sapply(cluster_id, function(id) {
  centers[id, ] + rnorm(2, sd = 0.3)
}))
y_cluster <- rnorm(n)

dcx_cluster <- build.nerve.from.knn(
  X = X_cluster,
  y = y_cluster,
  k = 20,
  max.p = 4,
  use.counting.measure = TRUE,
  directed.knn = FALSE
)

summary_cluster <- summary(dcx_cluster)
print(summary_cluster)

cat("\n=== Checking clustered data for duplicates ===\n")
for (p in 0:summary_cluster$max_dimension) {
  check_duplicates_manual(dcx_cluster, p)
}

# ==============================================================================
# Test 5: Extreme case - very high k
# ==============================================================================

cat("\n\n=== Test 5: High k value ===\n")
set.seed(99)
n <- 50
X_high_k <- matrix(rnorm(n * 2), n, 2)
y_high_k <- rnorm(n)

dcx_high_k <- build.nerve.from.knn(
  X = X_high_k,
  y = y_high_k,
  k = n - 5,  # Very high k
  max.p = 3,
  use.counting.measure = TRUE,
  directed.knn = FALSE
)

summary_high_k <- summary(dcx_high_k)
print(summary_high_k)

cat("\n=== Checking high-k case for duplicates ===\n")
for (p in 0:summary_high_k$max_dimension) {
  check_duplicates_manual(dcx_high_k, p)
}

cat("\n=== All manual tests complete ===\n")

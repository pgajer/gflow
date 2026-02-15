# Setup for interactive debugging sessions
# =============================================================================
# Usage:
#   - Set REBUILD = TRUE when you've modified C++ code
#   - Set REBUILD = FALSE to just load the installed package (faster)
# =============================================================================

REBUILD <- TRUE  # <-- Toggle this flag

if (REBUILD) {
  message("Rebuilding package from source...")
  devtools::clean_dll()
  devtools::document()
  devtools::load_all()
  message("✓ Package rebuilt and loaded")
} else {
  message("Loading installed package...")
  library(gflow)
  message("✓ Package loaded")
}

# =============================================================================
# Helper Functions
# =============================================================================

# Define the repeat operator FIRST (before it's used)
if (!exists("%R%")) {
  `%R%` <- function(x, n) paste(rep(x, n), collapse = "")
}

# Helper to reload after C++ changes (during same session)
reload_cpp <- function() {
  devtools::clean_dll()
  devtools::load_all()
  message("✓ Package reloaded")
}

# Quick test data generator
make_test_data <- function(n = 50, d = 2, seed = 123, type = "random") {
  set.seed(seed)

  X <- switch(type,
    random = matrix(rnorm(n * d), n, d),

    clustered = {
      # Create k clusters
      k_clusters <- ceiling(sqrt(n/10))
      centers <- matrix(rnorm(k_clusters * d, sd = 3), k_clusters, d)
      cluster_id <- sample(1:k_clusters, n, replace = TRUE)
      t(sapply(cluster_id, function(id) {
        centers[id, ] + rnorm(d, sd = 0.3)
      }))
    },

    grid = {
      # Regular grid with noise
      if (d != 2) stop("grid type only works for d=2")
      grid_size <- ceiling(sqrt(n))
      grid_points <- expand.grid(
        x = seq(0, 1, length.out = grid_size),
        y = seq(0, 1, length.out = grid_size)
      )
      as.matrix(grid_points[1:n, ]) + matrix(rnorm(n * 2, sd = 0.05), n, 2)
    },

    circle = {
      # Points on a circle with noise
      if (d != 2) stop("circle type only works for d=2")
      theta <- seq(0, 2*pi, length.out = n + 1)[1:n]
      cbind(cos(theta), sin(theta)) + matrix(rnorm(n * 2, sd = 0.1), n, 2)
    },

    stop("Unknown type. Use: random, clustered, grid, or circle")
  )

  # Generate response
  if (d == 2) {
    y <- X[, 1]^2 + X[, 2]^2 + rnorm(n, sd = 0.1)
  } else {
    y <- rowSums(X^2) + rnorm(n, sd = 0.1)
  }

  list(X = X, y = y, type = type, n = n, d = d)
}

# Print simplices with nice formatting
print_simplices <- function(dcx, dim, max_show = 10) {
  simplices <- .Call("S_get_simplices", dcx, as.integer(dim),
                     PACKAGE = "gflow")
  n <- length(simplices)
  cat(sprintf("\n=== Dimension %d: %d simplices ===\n", dim, n))

  if (n == 0) {
    cat("  (no simplices)\n")
    return(invisible(NULL))
  }

  show_n <- min(max_show, n)
  for (i in seq_len(show_n)) {
    cat(sprintf("  %4d: {%s}\n", i,
                paste(simplices[[i]], collapse = ", ")))
  }

  if (n > show_n) {
    cat(sprintf("  ... (%d more)\n", n - show_n))
  }

  invisible(simplices)
}

# Print summary of complex
print_complex_summary <- function(dcx) {
  summary_info <- riem.dcx.summary(dcx)

  cat("\n=== Nerve Complex Summary ===\n")
  cat(sprintf("Max dimension: %d\n", summary_info$max_dimension))
  cat(sprintf("Vertices:      %d\n", summary_info$n_vertices))
  cat(sprintf("Edges:         %d\n", summary_info$n_edges))
  if (summary_info$max_dimension >= 2) {
    cat(sprintf("Triangles:     %d\n", summary_info$n_triangles))
  }

  if (summary_info$max_dimension >= 3) {
    n_3 <- summary_info$n_simplices[4]  # 0-indexed in vector
    cat(sprintf("3-simplices:   %d\n", n_3))
  }

  cat(sprintf("Has response:  %s\n",
              ifelse(summary_info$has_response, "yes", "no")))

  invisible(summary_info)
}

# Check for duplicates across all dimensions
check_all_duplicates <- function(dcx, verbose = TRUE) {
  summary_info <- riem.dcx.summary(dcx)
  results <- list()

  if (verbose) cat("\n=== Checking for Duplicate Simplices ===\n")

  for (p in 0:summary_info$max_dimension) {
    simplices <- .Call("S_get_simplices", dcx, as.integer(p),
                       PACKAGE = "gflow")

    if (length(simplices) == 0) {
      if (verbose) cat(sprintf("Dimension %d: (empty)\n", p))
      results[[p + 1]] <- list(has_duplicates = FALSE, n_duplicates = 0)
      next
    }

    # Convert to strings for comparison
    simplex_strings <- vapply(simplices, function(s) {
      paste(sort(s), collapse = ",")
    }, character(1))

    # Find duplicates
    dup_mask <- duplicated(simplex_strings)
    has_dups <- any(dup_mask)
    n_dups <- sum(dup_mask)

    if (verbose) {
      if (has_dups) {
        cat(sprintf("Dimension %d: ✗ %d duplicates found!\n", p, n_dups))

        # Show first few
        dup_indices <- which(dup_mask)
        show_n <- min(3, n_dups)
        for (i in seq_len(show_n)) {
          idx <- dup_indices[i]
          cat(sprintf("  Duplicate %d: {%s}\n",
                      idx,
                      paste(simplices[[idx]], collapse = ", ")))
        }
        if (n_dups > show_n) {
          cat(sprintf("  ... (%d more)\n", n_dups - show_n))
        }
      } else {
        cat(sprintf("Dimension %d: ✓ No duplicates\n", p))
      }
    }

    results[[p + 1]] <- list(
      has_duplicates = has_dups,
      n_duplicates = n_dups,
      duplicate_indices = if (has_dups) which(dup_mask) else integer(0)
    )
  }

  invisible(results)
}

# Quick build and check
quick_check <- function(X = NULL, y = NULL, k = 10, max.p = 3, ...) {
  if (is.null(X)) {
    data <- make_test_data(n = 50, type = "random")
    X <- data$X
    y <- data$y
  }

  cat("\nBuilding nerve complex...\n")
  dcx <- build.nerve.from.knn(X = X, y = y, k = k, max.p = max.p, ...)

  print_complex_summary(dcx)
  check_all_duplicates(dcx)

  invisible(dcx)
}

# Visualize 2D complex (requires ggplot2)
plot_nerve_2d <- function(dcx, X, show_edges = TRUE, show_triangles = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting")
  }

  if (ncol(X) != 2) {
    stop("X must have exactly 2 columns for 2D plotting")
  }

  df <- data.frame(x = X[, 1], y = X[, 2], id = seq_len(nrow(X)))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Nerve Complex Visualization")

  # Add triangles (filled)
  if (show_triangles) {
    summary_info <- riem.dcx.summary(dcx)
    if (summary_info$max_dimension >= 2) {
      triangles <- .Call("S_get_simplices", dcx, 2L, PACKAGE = "gflow")

      if (length(triangles) > 0) {
        triangle_df <- do.call(rbind, lapply(triangles, function(tri) {
          data.frame(
            x = X[tri, 1],
            y = X[tri, 2],
            group = paste(tri, collapse = "-")
          )
        }))

        p <- p + ggplot2::geom_polygon(
          data = triangle_df,
          ggplot2::aes(x = x, y = y, group = group),
          fill = "lightblue",
          alpha = 0.3,
          color = NA
        )
      }
    }
  }

  # Add edges
  if (show_edges) {
    edges <- .Call("S_get_simplices", dcx, 1L, PACKAGE = "gflow")

    if (length(edges) > 0) {
      edge_df <- do.call(rbind, lapply(edges, function(e) {
        data.frame(
          x = X[e, 1],
          y = X[e, 2],
          group = paste(e, collapse = "-")
        )
      }))

      p <- p + ggplot2::geom_line(
        data = edge_df,
        ggplot2::aes(x = x, y = y, group = group),
        color = "gray50",
        alpha = 0.5
      )
    }
  }

  # Add vertices
  p <- p + ggplot2::geom_point(size = 2, color = "red")

  print(p)
  invisible(p)
}

# =============================================================================
# Print available functions
# =============================================================================

cat("\n")
cat("=" %R% 70, "\n")
cat("Debug Environment Ready\n")
cat("=" %R% 70, "\n")
cat("\nAvailable functions:\n")
cat("  Data generation:\n")
cat("    make_test_data(n, d, seed, type)  # type: random, clustered, grid, circle\n")
cat("\n")
cat("  Building and checking:\n")
cat("    quick_check(X, y, k, max.p, ...)  # Build and check for issues\n")
cat("    check_all_duplicates(dcx)         # Check all dimensions\n")
cat("\n")
cat("  Display functions:\n")
cat("    print_complex_summary(dcx)        # Overview of complex\n")
cat("    print_simplices(dcx, dim)         # Show simplices at dimension\n")
cat("    plot_nerve_2d(dcx, X)            # Visualize 2D complex (needs ggplot2)\n")
cat("\n")
cat("  Development:\n")
cat("    reload_cpp()                      # Reload after C++ changes\n")
cat("\n")
cat("Quick start:\n")
cat("  dcx <- quick_check()                # Run basic test\n")
cat("  data <- make_test_data(type='clustered'); dcx <- quick_check(data$X, data$y)\n")
cat("=" %R% 70, "\n\n")

test_that("build.nerve.from.knn constructs valid nerve complex", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  # Create test data
  set.seed(12345)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n))
  y <- X[, 1]^2 + X[, 2]^2 + rnorm(n, sd = 0.1)

  k <- 10
  max.p <- 3

  # Build nerve complex
  dcx <- build.nerve.from.knn(
    X = X,
    y = y,
    k = k,
    max.p = max.p,
    use.counting.measure = TRUE,
    directed.knn = FALSE
  )

  # Get summary
  summary.info <- riem.dcx.summary(dcx)

  # Basic structure checks
  expect_s3_class(dcx, "riem_dcx")
  expect_true(summary.info$max_dimension >= 1)
  expect_equal(summary.info$n_vertices, n)
  expect_true(summary.info$n_edges > 0)
})

test_that("nerve complex has no duplicate simplices", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  set.seed(12345)
  n <- 50
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

  summary_info <- riem.dcx.summary(dcx)
  max_dim <- summary_info$max_dimension

  # Check each dimension for duplicates
  for (p in 0:max_dim) {
    simplices <- get_simplices(dcx, p)

    if (length(simplices) > 0) {
      dup_check <- check_duplicates(simplices)

      expect_false(
        dup_check$has_duplicates,
        info = sprintf(
          "Dimension %d has duplicate simplices at indices: %s",
          p,
          paste(dup_check$duplicate_indices, collapse = ", ")
        )
      )

      # If duplicates found, show details for debugging
      if (dup_check$has_duplicates) {
        message(sprintf("Duplicate simplices in dimension %d:", p))
        for (i in seq_along(dup_check$duplicates)) {
          message(sprintf("  Simplex %d: {%s}",
                         dup_check$duplicate_indices[i],
                         paste(dup_check$duplicates[[i]], collapse = ", ")))
        }
      }
    }
  }
})

test_that("simplices are properly sorted", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  set.seed(12345)
  n <- 40
  X <- matrix(rnorm(n * 3), n, 3)
  y <- rnorm(n)

  dcx <- build.nerve.from.knn(X, y, k = 8, max.p = 2)
  summary_info <- riem.dcx.summary(dcx)

  for (p in 0:summary_info$max_dimension) {
    simplices <- get_simplices(dcx, p)

    if (length(simplices) > 0) {
      expect_true(
        simplices_are_sorted(simplices),
        info = sprintf("Simplices in dimension %d are not sorted", p)
      )
    }
  }
})

test_that("all simplex faces exist in lower dimension", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  set.seed(12345)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n))
  y <- rnorm(n)

  dcx <- build.nerve.from.knn(X, y, k = 10, max.p = 3)
  summary_info <- riem.dcx.summary(dcx)
  max_dim <- summary_info$max_dimension

  # Check face existence for dimensions >= 1
  for (p in 1:max_dim) {
    simplices_p <- get_simplices(dcx, p)
    simplices_p_minus_1 <- get_simplices(dcx, p - 1)

    if (length(simplices_p) > 0) {
      for (i in seq_along(simplices_p)) {
        has_all_faces <- all_faces_exist(simplices_p[[i]], simplices_p_minus_1)

        expect_true(
          has_all_faces,
          info = sprintf(
            "Simplex {%s} in dimension %d has missing faces",
            paste(simplices_p[[i]], collapse = ", "),
            p
          )
        )
      }
    }
  }
})

test_that("metric values are positive", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  set.seed(12345)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n))
  y <- rnorm(n)

  dcx <- build.nerve.from.knn(X, y, k = 10, max.p = 2)
  summary_info <- riem.dcx.summary(dcx)

  for (p in 0:summary_info$max_dimension) {
    metric_diag <- get_metric_diagonal(dcx, p)

    expect_true(
      all(metric_diag > 0),
      info = sprintf("Dimension %d has non-positive metric values", p)
    )

    expect_true(
      all(is.finite(metric_diag)),
      info = sprintf("Dimension %d has non-finite metric values", p)
    )
  }
})

test_that("directed vs undirected kNN give different results", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  set.seed(12345)
  n <- 40
  X <- cbind(rnorm(n), rnorm(n))
  y <- rnorm(n)

  dcx_directed <- build.nerve.from.knn(
    X, y, k = 10, max.p = 2,
    directed.knn = TRUE
  )

  dcx_undirected <- build.nerve.from.knn(
    X, y, k = 10, max.p = 2,
    directed.knn = FALSE
  )

  sum_dir <- riem.dcx.summary(dcx_directed)
  sum_undir <- riem.dcx.summary(dcx_undirected)

  # Directed should typically have more edges
  expect_true(
    sum_dir$n_edges >= sum_undir$n_edges,
    info = "Directed kNN should have at least as many edges as undirected (mutual) kNN"
  )
})

test_that("counting measure vs weighted measure behave correctly", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  set.seed(12345)
  n <- 40
  X <- cbind(rnorm(n), rnorm(n))
  y <- rnorm(n)

  dcx_counting <- build.nerve.from.knn(
    X, y, k = 10, max.p = 2,
    use.counting.measure = TRUE
  )

  dcx_weighted <- build.nerve.from.knn(
    X, y, k = 10, max.p = 2,
    use.counting.measure = FALSE
  )

  # Both should produce valid complexes
  expect_s3_class(dcx_counting, "riem_dcx")
  expect_s3_class(dcx_weighted, "riem_dcx")

  # Metric values should differ
  metric_count_0 <- get_metric_diagonal(dcx_counting, 0)
  metric_weight_0 <- get_metric_diagonal(dcx_weighted, 0)

  # They should not be identical (unless by extreme coincidence)
  expect_false(identical(metric_count_0, metric_weight_0))
})

test_that("complex construction handles edge cases", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  # Very small k
  set.seed(12345)
  n <- 20
  X <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)

  expect_no_error({
    dcx_small_k <- build.nerve.from.knn(X, y, k = 3, max.p = 2)
  })

  # k close to n
  expect_no_error({
    dcx_large_k <- build.nerve.from.knn(X, y, k = n - 2, max.p = 2)
  })
})

test_that("high dimensional simplices are constructed correctly", {
  skip_on_cran()
  skip_if_not(exists("build.nerve.from.knn", mode = "function") && exists("riem.dcx.summary", mode = "function"), "Legacy nerve API not available in this build")

  # Create data where we expect high-dimensional simplices
  set.seed(12345)
  n <- 100
  # Clustered data to increase overlap
  centers <- rbind(c(0, 0), c(3, 3), c(0, 3))
  cluster_id <- sample(1:3, n, replace = TRUE)
  X <- t(sapply(cluster_id, function(id) {
    centers[id, ] + rnorm(2, sd = 0.5)
  }))
  y <- rnorm(n)

  dcx <- build.nerve.from.knn(X, y, k = 15, max.p = 4)
  summary_info <- riem.dcx.summary(dcx)

  # Should have at least triangles
  expect_true(summary_info$max_dimension >= 2)

  if (summary_info$max_dimension >= 3) {
    # Check 3-simplices exist and have no duplicates
    simplices_3 <- get_simplices(dcx, 3)

    if (length(simplices_3) > 0) {
      dup_check <- check_duplicates(simplices_3)
      expect_false(dup_check$has_duplicates)

      # Each 3-simplex should have 4 vertices
      expect_true(all(vapply(simplices_3, length, integer(1)) == 4))
    }
  }
})

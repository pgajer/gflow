test_that("subject.neighborhood.stats computes expected diagnostics for conductance weights", {
  fit <- list(
    fitted.values = c(0.1, 0.2, 0.3, 0.4),
    graph = list(
      adj.list = list(
        c(2L, 3L),
        c(1L, 3L),
        c(1L, 2L, 4L),
        c(3L)
      ),
      edge.list = rbind(
        c(1L, 2L),
        c(1L, 3L),
        c(2L, 3L),
        c(3L, 4L)
      ),
      edge.densities = c(2, 4, 1, 0.5),
      vertex.densities = c(1, 4, 1, 2)
    )
  )

  subj.id <- c("A", "A", "B", "C")

  stats <- subject.neighborhood.stats(
    fit,
    subj.id,
    weight.type = "conductance",
    include.self.in.R = TRUE
  )

  expect_equal(stats$R, c(1.5, 1.5, 4 / 3, 1))
  expect_equal(stats$p.max, c(2 / 3, 2 / 3, 0.6153846154, 1), tolerance = 1e-10)
  expect_equal(stats$s.eff, c(1.8, 1.8, 1.8988764045, 1), tolerance = 1e-10)
})

test_that("subject.neighborhood.stats supports mass.sym weights and optimal.fit wrapper", {
  fit <- list(
    fitted.values = c(0.1, 0.2, 0.3, 0.4),
    graph = list(
      adj.list = list(
        c(2L, 3L),
        c(1L, 3L),
        c(1L, 2L, 4L),
        c(3L)
      ),
      edge.list = rbind(
        c(1L, 2L),
        c(1L, 3L),
        c(2L, 3L),
        c(3L, 4L)
      ),
      edge.densities = c(2, 4, 1, 0.5),
      vertex.densities = c(1, 4, 1, 2)
    )
  )

  wrapped <- list(optimal.fit = fit)
  subj.id <- c("A", "A", "B", "C")

  stats <- subject.neighborhood.stats(
    wrapped,
    subj.id,
    weight.type = "mass.sym",
    include.self.in.R = FALSE
  )

  expect_equal(stats$R, c(1, 1, 1.5, 1))
  expect_equal(stats$p.max, c(0.5, 2 / 3, 0.6534537935, 1), tolerance = 1e-8)
  expect_equal(stats$s.eff, c(2, 1.8, 1.8278312164, 1), tolerance = 1e-6)
})

test_that("subject.neighborhood.stats validates subject.id length", {
  fit <- list(
    fitted.values = c(0.1, 0.2),
    graph = list(
      adj.list = list(c(2L), c(1L)),
      edge.list = rbind(c(1L, 2L)),
      edge.densities = 1,
      vertex.densities = c(1, 1)
    )
  )

  expect_error(
    subject.neighborhood.stats(fit, subject.id = "s1"),
    "Length of subject.id"
  )
})

test_that("inverted-index graph builder matches pairscan reference", {
  set.seed(1)
  X <- matrix(rnorm(120 * 8), nrow = 120, ncol = 8)

  cmp <- .Call(
    "S_compare_iknn_graph_builders",
    X,
    as.integer(9L),
    as.integer(1L),
    as.logical(FALSE),
    PACKAGE = "gflow"
  )

  expect_true(is.list(cmp))
  expect_true(isTRUE(cmp$identical))
  expect_identical(cmp$n_missing_edges, 0L)
  expect_identical(cmp$n_extra_edges, 0L)
  expect_identical(cmp$n_isize_mismatch, 0L)
  expect_identical(cmp$n_dist_mismatch, 0L)
})

test_that("iknn graph constructors default to skipping optional stages", {
  set.seed(2)
  X <- matrix(rnorm(90 * 6), nrow = 90, ncol = 6)

  g.multi <- create.iknn.graphs(
    X,
    kmin = 4L,
    kmax = 5L,
    compute.full = TRUE,
    n.cores = 1L,
    verbose = FALSE
  )

  expect_null(g.multi$isize_pruned_graphs)
  expect_null(g.multi$edge_pruning_stats)
  expect_true(is.matrix(g.multi$k_statistics))
  expect_true(all(is.na(g.multi$k_statistics[, "n_edges_in_isize_pruned_graph"])))

  g.single <- create.single.iknn.graph(
    X,
    k = 5L,
    compute.full = FALSE,
    verbose = FALSE
  )

  expect_null(g.single$isize_pruned_adj_list)
  expect_null(g.single$isize_pruned_weight_list)
  expect_true(is.na(g.single$n_edges_in_isize_pruned_graph))
  expect_null(g.single$edge_pruning_stats)
})

test_that("create.single.iknn.graph supports strict kNN cache read/write", {
  set.seed(3)
  X <- matrix(rnorm(80 * 5), nrow = 80, ncol = 5)
  cache_path <- tempfile("iknn_single_cache_", fileext = ".bin")
  on.exit(unlink(cache_path), add = TRUE)

  g_write <- create.single.iknn.graph(
    X,
    k = 6L,
    compute.full = FALSE,
    knn.cache.path = cache_path,
    knn.cache.mode = "write",
    verbose = FALSE
  )
  expect_true(file.exists(cache_path))

  g_read <- create.single.iknn.graph(
    X,
    k = 6L,
    compute.full = FALSE,
    knn.cache.path = cache_path,
    knn.cache.mode = "read",
    verbose = FALSE
  )

  expect_equal(g_write$pruned_adj_list, g_read$pruned_adj_list)
  expect_equal(g_write$pruned_weight_list, g_read$pruned_weight_list, tolerance = 1e-12)
})

test_that("kNN cache validation rejects dimension mismatches", {
  set.seed(4)
  X <- matrix(rnorm(70 * 4), nrow = 70, ncol = 4)
  X_bad <- cbind(X, rnorm(nrow(X)))
  cache_path <- tempfile("iknn_dim_mismatch_", fileext = ".bin")
  on.exit(unlink(cache_path), add = TRUE)

  invisible(create.single.iknn.graph(
    X,
    k = 5L,
    compute.full = FALSE,
    knn.cache.path = cache_path,
    knn.cache.mode = "write",
    verbose = FALSE
  ))

  expect_error(
    create.single.iknn.graph(
      X_bad,
      k = 5L,
      compute.full = FALSE,
      knn.cache.path = cache_path,
      knn.cache.mode = "read",
      verbose = FALSE
    ),
    "ncol mismatch"
  )
})

test_that("create.iknn.graphs readwrite cache computes on miss and saves", {
  set.seed(5)
  X <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)
  cache_path <- tempfile("iknn_multi_cache_", fileext = ".bin")
  on.exit(unlink(cache_path), add = TRUE)

  expect_false(file.exists(cache_path))

  g <- create.iknn.graphs(
    X,
    kmin = 4L,
    kmax = 6L,
    compute.full = FALSE,
    n.cores = 1L,
    knn.cache.path = cache_path,
    knn.cache.mode = "readwrite",
    verbose = FALSE
  )

  expect_true(file.exists(cache_path))
  expect_true(is.matrix(g$k_statistics))
  expect_equal(nrow(g$k_statistics), 3L)
})

test_that("kNN cache write expands tilde and creates parent directories", {
  set.seed(6)
  X <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
  fake_home <- tempfile("iknn_fake_home_")
  old_home <- Sys.getenv("HOME", unset = "")
  on.exit(Sys.setenv(HOME = old_home), add = TRUE)
  on.exit(unlink(fake_home, recursive = TRUE), add = TRUE)
  dir.create(fake_home, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(HOME = fake_home)

  cache_path <- "~/nested/cache/knn.bin"
  expanded_cache_path <- path.expand(cache_path)
  expect_false(file.exists(expanded_cache_path))

  invisible(create.single.iknn.graph(
    X,
    k = 6L,
    compute.full = FALSE,
    knn.cache.path = cache_path,
    knn.cache.mode = "write",
    verbose = FALSE
  ))

  expect_true(file.exists(expanded_cache_path))
})

test_that("kNN cache path rejects directory input", {
  set.seed(7)
  X <- matrix(rnorm(40 * 4), nrow = 40, ncol = 4)
  cache_dir <- tempfile("iknn_cache_dir_")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  expect_error(
    create.single.iknn.graph(
      X,
      k = 5L,
      compute.full = FALSE,
      knn.cache.path = cache_dir,
      knn.cache.mode = "write",
      verbose = FALSE
    ),
    "cache filename"
  )
})

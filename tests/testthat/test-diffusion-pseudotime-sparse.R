make_path_graph <- function(n = 120L, k = 4L) {
  stopifnot(n > 10, k >= 1)
  adj <- vector("list", n)
  w <- vector("list", n)
  x <- seq_len(n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - k)
    hi <- min(n, i + k)
    nbr <- setdiff(lo:hi, i)
    dist <- abs(x[nbr] - x[i])
    adj[[i]] <- as.integer(nbr)
    w[[i]] <- as.numeric(dist + 1)
  }
  list(adj = adj, w = w)
}


test_that("potential sparse pseudotime matches dense matrix-power reference", {
  skip_if_not_installed("Matrix")

  g <- make_path_graph(n = 140L, k = 3L)
  roots <- c(1L, 2L, 3L, 4L)
  t.steps <- 5L
  eps <- 1e-10

  tr <- build.sparse.transition(
    adj.list = g$adj,
    weight.list = g$w,
    weight.mode = "inverse",
    weight.param = 1e-6,
    lazy = 1
  )
  P <- as.matrix(tr$P)
  Pt <- P
  for (ii in 2:t.steps) {
    Pt <- Pt %*% P
  }

  root.mix <- rep(0, nrow(P))
  root.mix[roots] <- 1 / length(roots)
  score.ref <- as.numeric(t(Pt) %*% root.mix)
  pot.ref <- -log(pmax(score.ref, eps))
  pot.ref.norm <- (pot.ref - min(pot.ref)) / (max(pot.ref) - min(pot.ref))

  fit <- compute.potential.pseudotime.sparse(
    adj.list = g$adj,
    weight.list = g$w,
    root.vertices = roots,
    t.steps = t.steps,
    potential.eps = eps,
    weight.mode = "inverse",
    weight.param = 1e-6,
    normalize = TRUE,
    return.transition = TRUE
  )

  expect_s3_class(fit, "potential_pseudotime_sparse")
  expect_equal(length(fit$pseudotime), length(score.ref))
  expect_equal(length(fit$root.score), length(score.ref))
  expect_equal(length(fit$root.potential), length(score.ref))
  expect_s4_class(fit$transition, "dgCMatrix")

  expect_equal(fit$root.score, score.ref, tolerance = 1e-10)
  expect_equal(fit$root.potential, pot.ref, tolerance = 1e-10)
  expect_equal(fit$pseudotime, pot.ref.norm, tolerance = 1e-10)
})


test_that("diffusion sparse pseudotime matches dense rooted diffusion-score baseline", {
  g <- make_path_graph(n = 120L, k = 3L)
  roots <- c(1L, 2L, 3L, 4L)
  t.steps <- 4L

  tr <- build.sparse.transition(
    adj.list = g$adj,
    weight.list = g$w,
    weight.mode = "inverse",
    weight.param = 1e-6,
    lazy = 1
  )
  P <- as.matrix(tr$P)
  Pt <- P
  for (ii in 2:t.steps) {
    Pt <- Pt %*% P
  }
  root.mix <- rep(0, nrow(P))
  root.mix[roots] <- 1 / length(roots)
  score.ref <- as.numeric(t(Pt) %*% root.mix)
  dist.ref <- max(score.ref) - score.ref

  fit <- compute.diffusion.pseudotime.sparse(
    adj.list = g$adj,
    weight.list = g$w,
    root.vertices = roots,
    t.steps = t.steps,
    n.probes = 256L,
    seed = 7L,
    weight.mode = "inverse",
    weight.param = 1e-6,
    normalize = FALSE
  )

  expect_s3_class(fit, "diffusion_pseudotime_sparse")
  expect_equal(length(fit$diffusion.distance), length(dist.ref))
  expect_equal(length(fit$pseudotime), length(dist.ref))
  expect_equal(length(fit$root.score), length(score.ref))

  expect_equal(fit$root.score, score.ref, tolerance = 1e-10)
  expect_equal(fit$diffusion.distance, dist.ref, tolerance = 1e-10)
})


test_that("potential sparse landmark mode returns embedding and optional distances", {
  g <- make_path_graph(n = 90L, k = 2L)
  roots <- c(1L, 2L, 3L)
  landmarks <- c(1L, 15L, 30L, 45L, 60L, 75L, 90L)

  fit <- compute.potential.pseudotime.sparse(
    adj.list = g$adj,
    weight.list = g$w,
    root.vertices = roots,
    t.steps = 4L,
    landmark.vertices = landmarks,
    return.landmark.distances = TRUE,
    max.landmark.distance.n = 200L,
    weight.mode = "inverse",
    weight.param = 1e-6,
    normalize = TRUE
  )

  expect_true(!is.null(fit$landmark.potential))
  expect_equal(dim(fit$landmark.potential), c(length(g$adj), length(landmarks)))
  expect_true(!is.null(fit$landmark.distances))
  expect_equal(dim(fit$landmark.distances), c(length(g$adj), length(g$adj)))
})


test_that("potential sparse root score matches phate.core dense Pt aggregation", {
  g <- make_path_graph(n = 100L, k = 3L)
  roots <- c(1L, 2L, 3L, 4L, 5L)
  t.steps <- 4L
  eps <- 1e-12

  tr <- build.sparse.transition(
    adj.list = g$adj,
    weight.list = g$w,
    weight.mode = "inverse",
    weight.param = 1e-6,
    lazy = 1
  )
  P <- as.matrix(tr$P)

  core <- phate.core(
    P = P,
    t = t.steps,
    potential.eps = eps,
    compute.D.pot = FALSE,
    verbose = FALSE
  )
  score.ref <- colMeans(core$Pt[roots, , drop = FALSE])

  fit <- compute.potential.pseudotime.sparse(
    adj.list = g$adj,
    weight.list = g$w,
    root.vertices = roots,
    t.steps = t.steps,
    potential.eps = eps,
    weight.mode = "inverse",
    weight.param = 1e-6,
    normalize = FALSE
  )

  expect_equal(fit$root.score, as.numeric(score.ref), tolerance = 1e-10)
})

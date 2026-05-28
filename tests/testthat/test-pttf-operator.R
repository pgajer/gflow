make_pttf_operator_path_graph <- function(n) {
  adj <- vector("list", n)
  w <- vector("list", n)
  x <- seq(0, 1, length.out = n)
  for (i in seq_len(n - 1L)) {
    ell <- x[i + 1L] - x[i]
    adj[[i]] <- c(adj[[i]], i + 1L)
    w[[i]] <- c(w[[i]], ell)
    adj[[i + 1L]] <- c(adj[[i + 1L]], i)
    w[[i + 1L]] <- c(w[[i + 1L]], ell)
  }
  list(
    X = matrix(x, ncol = 1L),
    adj.list = lapply(adj, as.integer),
    weight.list = lapply(w, as.double)
  )
}

make_pttf_operator_path_geometry <- function(n = 12L) {
  graph <- make_pttf_operator_path_graph(n)
  geom <- pttf.geometry(
    graph$X,
    graph = "supplied",
    adj.list = graph$adj.list,
    weight.list = graph$weight.list,
    tangent.dim = 1L
  )
  list(graph = graph, geom = geom)
}

test_that("PTTF symmetric tensor bases have expected dimensions", {
  expect_equal(gflow:::.pttf.operator.q(1L, 3L), 1L)
  expect_equal(gflow:::.pttf.operator.q(2L, 3L), 4L)
  expect_equal(gflow:::.pttf.operator.q(3L, 3L), 10L)

  b2 <- gflow:::.pttf.symmetric.tensor.basis(3L, 2L, "hs")
  expect_equal(b2$q, 6L)
  expect_equal(b2$multiplicity, c(1L, 2L, 2L, 1L, 2L, 1L))

  b3 <- gflow:::.pttf.symmetric.tensor.basis(3L, 3L, "hs")
  expect_equal(b3$q, 10L)
  expect_true(all(b3$multiplicity %in% c(1L, 3L, 6L)))
  expect_equal(sum(b3$multiplicity), 27L)
})

test_that("PTTF tensor design contracts symmetric tensors correctly", {
  u <- c(0.5, -1.25)
  H <- matrix(c(2, 3, 3, -4), nrow = 2L)
  T <- array(0, dim = c(2L, 2L, 2L))
  T[1, 1, 1] <- 1
  T[1, 1, 2] <- T[1, 2, 1] <- T[2, 1, 1] <- 2
  T[1, 2, 2] <- T[2, 1, 2] <- T[2, 2, 1] <- -1
  T[2, 2, 2] <- 0.5

  for (scaling in c("hs", "raw")) {
    b1 <- gflow:::.pttf.symmetric.tensor.basis(2L, 1L, scaling)
    b2 <- gflow:::.pttf.symmetric.tensor.basis(2L, 2L, scaling)
    b3 <- gflow:::.pttf.symmetric.tensor.basis(2L, 3L, scaling)

    h <- gflow:::.pttf.tensor.coord.from.full(H, b2)
    phi2 <- gflow:::.pttf.symmetric.tensor.direction.design(u, 1L, 2L, scaling)
    expect_equal(as.vector(phi2 %*% h),
                 gflow:::.pttf.tensor.coord.from.full(as.vector(H %*% u), b1),
                 tolerance = 1e-12,
                 info = scaling)

    tcoord <- gflow:::.pttf.tensor.coord.from.full(T, b3)
    phi3 <- gflow:::.pttf.symmetric.tensor.direction.design(u, 2L, 3L, scaling)
    contracted <- gflow:::.pttf.tensor.contract.last(T, u)
    expect_equal(as.vector(phi3 %*% tcoord),
                 gflow:::.pttf.tensor.coord.from.full(contracted, b2),
                 tolerance = 1e-12,
                 info = scaling)
  }
})

test_that("PTTF tensor transport matches explicit matrix formulas", {
  angle <- 0.37
  O <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2L)
  H <- matrix(c(1.5, -0.4, -0.4, 2.5), nrow = 2L)

  for (scaling in c("hs", "raw")) {
    b2 <- gflow:::.pttf.symmetric.tensor.basis(2L, 2L, scaling)
    h <- gflow:::.pttf.tensor.coord.from.full(H, b2)
    K2 <- gflow:::.pttf.symmetric.tensor.transport(O, 2L, scaling)
    expect_equal(as.vector(K2 %*% h),
                 gflow:::.pttf.tensor.coord.from.full(O %*% H %*% t(O), b2),
                 tolerance = 1e-12,
                 info = scaling)
  }
})

test_that("pttf.operator assembles compact sparse operators and provenance", {
  obj <- make_pttf_operator_path_geometry(12L)
  op1 <- pttf.operator(obj$geom, derivative.order = 1L,
                       row.mass.rule = "none", row.normalize = "none")
  expect_s3_class(op1, "pttf_operator")
  expect_s4_class(op1$A, "dgCMatrix")
  expect_equal(dim(op1$A), c(12L, 12L))
  expect_equal(dim(op1$B), c(12L, 12L))
  expect_equal(as.matrix(op1$B), as.matrix(Matrix::crossprod(op1$A)),
               tolerance = 1e-12)
  expect_equal(nrow(op1$row.table), nrow(op1$A))
  expect_true(all(op1$row.table$full.row == seq_len(nrow(op1$A))))
  expect_true(all(op1$vertex.table$edge.policy == "ok.only"))
  expect_equal(max(abs(as.vector(op1$A %*% rep(1, 12L)))), 0,
               tolerance = 1e-10)
})

test_that("pttf.operator tracks predecessor availability instead of zero filling", {
  obj <- make_pttf_operator_path_geometry(12L)
  geom <- obj$geom
  geom$frames$support.diagnostics$status[6L] <- "rank_deficient"
  op2 <- pttf.operator(geom, derivative.order = 2L,
                       row.mass.rule = "none", row.normalize = "none")
  level2 <- op2$vertex.table[op2$vertex.table$derivative.order == 2L, ]
  expect_true(any(level2$status == "predecessor_block_unavailable"))
  expect_true(any(level2$n.predecessor.neighbors.dropped > 0L))
  expect_true(all(op2$diagnostics$availability$available.level1[6L] == FALSE))
})

test_that("pttf.operator order-2 path ladder annihilates affine probes", {
  obj <- make_pttf_operator_path_geometry(12L)
  op2 <- pttf.operator(obj$geom, derivative.order = 2L,
                       row.mass.rule = "none", row.normalize = "none")
  x <- obj$graph$X[, 1L]
  probes <- cbind(1, x)
  residual <- as.matrix(op2$A %*% probes)
  expect_lt(max(abs(residual)), 1e-8)
})

test_that("pttf.operator order-3 path ladder annihilates quadratic probes away from boundaries", {
  obj <- make_pttf_operator_path_geometry(24L)
  op3 <- pttf.operator(obj$geom, derivative.order = 3L,
                       row.mass.rule = "none", row.normalize = "none")
  x <- obj$graph$X[, 1L]
  residual <- as.vector(op3$A %*% (x^2))
  interior <- op3$row.table$vertex > 3L & op3$row.table$vertex <= length(x) - 3L
  boundary <- !interior

  expect_lt(max(abs(residual[interior])), 1e-8)
  expect_gt(max(abs(residual[boundary])), 1e-4)
})

test_that("pttf.operator validates edge policy and row normalization options", {
  obj <- make_pttf_operator_path_geometry(12L)
  geom <- obj$geom
  geom$edges$status[5L] <- "too_few_shared_points"

  op.ok <- pttf.operator(geom, derivative.order = 1L,
                         edge.status.policy = "ok.only",
                         row.mass.rule = "none", row.normalize = "none")
  op.fallback <- pttf.operator(geom, derivative.order = 1L,
                               edge.status.policy = "frame.fallback",
                               row.mass.rule = "none", row.normalize = "l2")
  expect_gt(sum(op.fallback$vertex.table$usable.degree),
            sum(op.ok$vertex.table$usable.degree))
  expect_equal(max(abs(op.fallback$row.table$row.norm.after.normalization - 1)),
               0, tolerance = 1e-10)
})

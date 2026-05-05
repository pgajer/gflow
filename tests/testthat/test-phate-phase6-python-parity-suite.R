phase6_fixture_root <- function() {
  testthat::test_path("fixtures", "ph6", "core")
}

phase6_case_ids <- function() {
  root <- phase6_fixture_root()
  basename(list.dirs(root, recursive = FALSE, full.names = TRUE))
}

phase6_read_matrix <- function(path) {
  unname(as.matrix(utils::read.csv(path, header = FALSE)))
}

phase6_read_metadata_number <- function(case.dir, key) {
  lines <- readLines(file.path(case.dir, "metadata.json"), warn = FALSE)
  line <- grep(sprintf('"%s"', key), lines, value = TRUE)[1]
  as.numeric(sub(sprintf('.*"%s": ([0-9.]+).*', key), "\\1", line))
}

phase6_embedding_file <- function(case.dir, method) {
  suffix <- switch(method,
    classic = "c",
    metric = "m",
    nonmetric = "nm",
    stop("unknown embedding method")
  )
  file.path(case.dir, sprintf("e_%s.csv", suffix))
}

phase6_pairwise <- function(X) {
  as.numeric(stats::dist(X))
}

phase6_embedding_metrics <- function(ref, obs) {
  d.ref <- phase6_pairwise(ref)
  d.obs <- phase6_pairwise(obs)
  scale.factor <- sum(d.ref * d.obs) / sum(d.obs^2)
  rel.dist.error <- sqrt(sum((d.ref - scale.factor * d.obs)^2) / sum(d.ref^2))

  A <- scale(ref, center = TRUE, scale = FALSE)
  B <- scale(obs, center = TRUE, scale = FALSE)
  sv <- svd(t(B) %*% A)
  rotation <- sv$u %*% t(sv$v)
  B.rot <- B %*% rotation
  B.scale <- sum(A * B.rot) / sum(B.rot^2)
  procrustes.error <- sqrt(sum((A - B.scale * B.rot)^2) / sum(A^2))

  list(
    rel_dist_error = rel.dist.error,
    procrustes_error = procrustes.error,
    pearson = stats::cor(d.ref, d.obs, method = "pearson"),
    spearman = stats::cor(d.ref, d.obs, method = "spearman")
  )
}

test_that("PHATE Phase 6 core operators match frozen original-implementation fixtures", {
  for (case.id in phase6_case_ids()) {
    case.dir <- file.path(phase6_fixture_root(), case.id)
    X <- phase6_read_matrix(file.path(case.dir, "X.csv"))
    t.auto <- phase6_read_metadata_number(case.dir, "t_auto")

    fit <- phate.core(
      X = X,
      k = 5,
      alpha.decay = 40,
      t = "auto",
      t.max = 30,
      gamma = 1,
      kernel.mode = "phate",
      vne.method = "phate",
      potential.mode = "phate",
      compute.D.pot = TRUE,
      verbose = FALSE
    )

    K.ref <- phase6_read_matrix(file.path(case.dir, "K.csv"))
    P.ref <- phase6_read_matrix(file.path(case.dir, "P.csv"))
    Pt.ref <- phase6_read_matrix(file.path(case.dir, "Pt.csv"))
    U.ref <- phase6_read_matrix(file.path(case.dir, "U.csv"))
    D.ref <- phase6_read_matrix(file.path(case.dir, "D.csv"))
    vne.ref <- as.numeric(phase6_read_matrix(file.path(case.dir, "vne.csv")))

    expect_equal(fit$t, t.auto, info = case.id)
    expect_equal(fit$K, K.ref, tolerance = 1e-10, info = case.id)
    expect_equal(fit$P, P.ref, tolerance = 1e-10, info = case.id)
    expect_equal(fit$diagnostics$vne, vne.ref, tolerance = 1e-10, info = case.id)
    expect_equal(fit$Pt, Pt.ref, tolerance = 1e-10, info = case.id)
    expect_equal(fit$U.pot, U.ref, tolerance = 1e-10, info = case.id)
    expect_equal(unname(fit$D.pot), D.ref, tolerance = 1e-10, info = case.id)
  }
})

test_that("PHATE Phase 6 core embeddings match frozen fixtures geometrically", {
  for (case.id in phase6_case_ids()) {
    case.dir <- file.path(phase6_fixture_root(), case.id)
    X <- phase6_read_matrix(file.path(case.dir, "X.csv"))
    seed <- phase6_read_metadata_number(case.dir, "seed")

    core <- phate.core(
      X = X,
      k = 5,
      alpha.decay = 40,
      t = "auto",
      t.max = 30,
      gamma = 1,
      kernel.mode = "phate",
      vne.method = "phate",
      potential.mode = "phate",
      compute.D.pot = TRUE,
      verbose = FALSE
    )

    for (method in c("classic", "metric", "nonmetric")) {
      emb <- phate.embed(
        core = core,
        ndim = 2,
        method = method,
        seed = seed,
        maxit = 300,
        verbose = FALSE
      )
      ref <- phase6_read_matrix(phase6_embedding_file(case.dir, method))
      metrics <- phase6_embedding_metrics(ref, emb$embedding)
      label <- sprintf("%s/%s", case.id, method)

      expect_true(metrics$rel_dist_error <= 0.18, info = label)
      expect_true(metrics$procrustes_error <= 0.28, info = label)
      expect_true(metrics$pearson >= 0.94, info = label)
      expect_true(metrics$spearman >= 0.90, info = label)
    }
  }
})

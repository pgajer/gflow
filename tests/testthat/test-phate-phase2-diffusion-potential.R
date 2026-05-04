test_that("PHATE Phase 2 phate VNE auto-t matches reference", {
  phase1.dir <- testthat::test_path("fixtures", "phate_phase1")
  phase2.dir <- testthat::test_path("fixtures", "phate_phase2")

  P.ref <- unname(as.matrix(utils::read.csv(file.path(phase1.dir, "dense_36x5_P_python.csv"),
                                            header = FALSE)))
  vne.ref <- as.numeric(utils::read.csv(file.path(phase2.dir, "dense_36x5_vne_tmax30_python.csv"),
                                        header = FALSE)[1, ])
  meta <- readLines(file.path(phase2.dir, "dense_36x5_phase2_meta.txt"),
                    warn = FALSE)
  t.auto <- as.integer(sub("^t_auto=", "", meta[grep("^t_auto=", meta)]))

  fit <- phate.core(
    P = P.ref,
    t = "auto",
    t.max = 30,
    vne.method = "phate",
    gamma = 1,
    potential.mode = "phate",
    compute.D.pot = FALSE,
    verbose = FALSE
  )

  expect_equal(fit$diagnostics$vne_method, "phate")
  expect_equal(fit$diagnostics$t_grid, 0:29)
  expect_equal(fit$diagnostics$vne, vne.ref, tolerance = 1e-12)
  expect_equal(fit$t, t.auto)
})


test_that("PHATE Phase 2 phate diffusion potentials match gamma references", {
  phase1.dir <- testthat::test_path("fixtures", "phate_phase1")
  phase2.dir <- testthat::test_path("fixtures", "phate_phase2")

  P.ref <- unname(as.matrix(utils::read.csv(file.path(phase1.dir, "dense_36x5_P_python.csv"),
                                            header = FALSE)))
  Pt.ref <- unname(as.matrix(utils::read.csv(file.path(phase2.dir, "dense_36x5_Pt_t3_python.csv"),
                                             header = FALSE)))

  cases <- list(
    list(gamma = 1, suffix = "gamma1"),
    list(gamma = 0, suffix = "gamma0"),
    list(gamma = -1, suffix = "gamma_neg1")
  )

  for (case in cases) {
    U.ref <- unname(as.matrix(utils::read.csv(
      file.path(phase2.dir, sprintf("dense_36x5_U_t3_%s_python.csv", case$suffix)),
      header = FALSE
    )))
    D.ref <- unname(as.matrix(utils::read.csv(
      file.path(phase2.dir, sprintf("dense_36x5_Dpot_t3_%s_python.csv", case$suffix)),
      header = FALSE
    )))

    fit <- phate.core(
      P = P.ref,
      t = 3,
      gamma = case$gamma,
      potential.mode = "phate",
      compute.D.pot = TRUE,
      verbose = FALSE
    )

    expect_equal(fit$diagnostics$potential_mode, "phate")
    expect_equal(fit$diagnostics$potential_eps, 1e-7)
    expect_equal(fit$diagnostics$gamma, case$gamma)
    expect_equal(fit$Pt, Pt.ref, tolerance = 1e-12)
    expect_equal(fit$U.pot, U.ref, tolerance = 1e-10)
    expect_equal(unname(fit$D.pot), D.ref, tolerance = 1e-10)
  }
})


test_that("PHATE Phase 2 deprecated python aliases map to phate", {
  phase1.dir <- testthat::test_path("fixtures", "phate_phase1")
  P.ref <- unname(as.matrix(utils::read.csv(file.path(phase1.dir, "dense_36x5_P_python.csv"),
                                            header = FALSE)))

  expect_warning(
    fit.vne <- phate.core(
      P = P.ref,
      t = "auto",
      t.max = 30,
      vne.method = "python",
      gamma = 1,
      potential.mode = "phate",
      compute.D.pot = FALSE,
      verbose = FALSE
    ),
    "deprecated"
  )

  expect_warning(
    fit.potential <- phate.core(
      P = P.ref,
      t = 3,
      vne.method = "phate",
      gamma = 1,
      potential.mode = "python",
      compute.D.pot = FALSE,
      verbose = FALSE
    ),
    "deprecated"
  )

  expect_equal(fit.vne$diagnostics$vne_method, "phate")
  expect_equal(fit.potential$diagnostics$potential_mode, "phate")
})


test_that("PHATE Phase 2 preserves legacy gflow log-potential default", {
  P <- matrix(c(
    0.8, 0.2,
    0.0, 1.0
  ), nrow = 2, byrow = TRUE)

  fit <- phate.core(
    P = P,
    t = 1,
    gamma = 1,
    compute.D.pot = FALSE,
    verbose = FALSE
  )

  expect_equal(fit$diagnostics$potential_mode, "gflow")
  expect_equal(fit$diagnostics$potential_eps, 1e-12)
  expect_equal(fit$U.pot, -log(pmax(P, 1e-12)), tolerance = 1e-12)
})

posterior.path.fixture <- function() {
    adjacency <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
    list(
        adjacency = adjacency,
        lengths = lapply(adjacency, function(v) rep(1, length(v))),
        field = c(0, 3, 1, 2, 0)
    )
}

test_that("every feature is validated before any local correlation is computed", {
    f <- posterior.path.fixture()
    draws <- cbind(f$field, f$field)
    # An invalid graph would fail first if any feature were computed before
    # the later feature's input was validated.
    for (later in list(cbind(draws, -f$field), draws[, 1, drop = FALSE])) {
        expect_error(
            lcor.with.posterior(NULL, NULL, f$field,
                                list(aligned = draws, mixed = later),
                                verbose = FALSE),
            "Z.hat.samples\\[\\[2\\]\\].*mixed.*columns.*expected 2"
        )
    }
    invalid <- list(
        vector = f$field,
        rows = draws[-1, , drop = FALSE],
        character = matrix("bad", 5, 2),
        complex = draws + 1i,
        missing = replace(draws, 1, NA_real_),
        infinite = replace(draws, 1, Inf)
    )
    for (later in invalid) {
        expect_error(
            lcor.with.posterior(NULL, NULL, f$field,
                                list(aligned = draws, later = later),
                                verbose = FALSE),
            "Z.hat.samples\\[\\[2\\]\\].*later"
        )
    }
})

test_that("posterior summaries use all draws and match an independent formula", {
    f <- posterior.path.fixture()
    draws <- cbind(f$field, f$field, -f$field)
    nonlinear <- cbind(c(1, 0, 2, 4, 3), c(0, 2, 3, 1, 4), -f$field)
    inputs <- list(draws, nonlinear)
    result <- lcor.with.posterior(f$adjacency, f$lengths, f$field, inputs,
                                  credible.level = 0.8, return.samples = TRUE,
                                  verbose = FALSE)
    expect_equal(as.numeric(result$mean[1, ]), rep(1 / 3, 5))
    expect_identical(result$n.samples, 3L)
    for (j in seq_along(inputs)) {
        expected <- vapply(seq_len(ncol(inputs[[j]])), function(b) {
            vapply(seq_along(f$field), function(i) {
                dy <- f$field[f$adjacency[[i]]] - f$field[i]
                dz <- inputs[[j]][f$adjacency[[i]], b] - inputs[[j]][i, b]
                weights <- 1 / f$lengths[[i]]^2
                sum(weights * dy * dz) /
                    sqrt(sum(weights * dy^2) * sum(weights * dz^2))
            }, numeric(1))
        }, numeric(5))
        expect_equal(result$samples[[j]], expected, tolerance = 1e-12)
        expect_equal(as.numeric(result$mean[j, ]), rowMeans(expected))
        expect_equal(as.numeric(result$sd[j, ]), apply(expected, 1, sd))
        expect_equal(as.numeric(result$lower[j, ]),
                     apply(expected, 1, quantile, probs = 0.1))
        expect_equal(as.numeric(result$upper[j, ]),
                     apply(expected, 1, quantile, probs = 0.9))
    }
    single <- lcor.with.posterior(f$adjacency, f$lengths, f$field, draws,
                                  verbose = FALSE)
    expect_equal(single$mean, result$mean[1, , drop = FALSE])
    expect_null(dimnames(single$mean))
    expect_null(single$samples)
})

test_that("empty inputs, nonfinite fields, and invalid controls are rejected", {
    f <- posterior.path.fixture()
    draws <- cbind(f$field, -f$field)
    run <- function(z = draws, y = f$field, ...) {
        lcor.with.posterior(f$adjacency, f$lengths, y, z, ...)
    }
    for (z in list(list(), NULL, f$field, as.data.frame(draws), matrix(numeric(), 5, 0))) {
        expect_error(run(z, verbose = FALSE), "Z.hat.samples")
    }
    for (y in list(numeric(), as.character(f$field), matrix(f$field),
                   replace(f$field, 1, NA_real_), replace(f$field, 1, Inf))) {
        expect_error(run(y = y, verbose = FALSE), "y.hat")
    }
    for (level in list(0, 1, -0.1, 1.1, NA_real_, NaN, Inf, NULL, c(0.8, 0.9), "0.95")) {
        expect_error(run(credible.level = level, verbose = FALSE), "credible.level")
    }
    for (flag in list(NA, NULL, c(TRUE, FALSE), 1)) {
        expect_error(run(return.samples = flag, verbose = FALSE), "return.samples")
        expect_error(run(verbose = flag), "verbose")
    }
})

test_that("one draw has its observed coefficient as interval and undefined SD", {
    f <- posterior.path.fixture()
    result <- lcor.with.posterior(f$adjacency, f$lengths, f$field,
                                  matrix(f$field, ncol = 1), verbose = FALSE)
    expect_equal(as.numeric(result$mean), rep(1, 5))
    expect_identical(result$lower, result$mean)
    expect_identical(result$upper, result$mean)
    expect_true(all(is.na(result$sd)))
    expect_identical(result$n.samples, 1L)
})

test_that("feature and vertex identities survive summaries and saved samples", {
    f <- posterior.path.fixture()
    ids <- paste0("sample-", c(5, 3, 1, 4, 2))
    draws <- cbind(positive = f$field, negative = -f$field)
    rownames(draws) <- ids
    result <- lcor.with.posterior(f$adjacency, f$lengths,
                                  setNames(f$field, ids),
                                  list(marker = draws, other = draws),
                                  return.samples = TRUE, verbose = FALSE)
    for (field in c("mean", "sd", "lower", "upper")) {
        expect_identical(dimnames(result[[field]]), list(c("marker", "other"), ids))
    }
    expect_identical(names(result$samples), c("marker", "other"))
    expect_identical(dimnames(result$samples$marker), dimnames(draws))
    expect_setequal(summary(result)$feature.summary$feature, c("marker", "other"))
    from.rows <- lcor.with.posterior(f$adjacency, f$lengths, f$field, draws,
                                     verbose = FALSE)
    expect_identical(colnames(from.rows$mean), ids)
    from.response <- lcor.with.posterior(f$adjacency, f$lengths,
                                         setNames(f$field, ids), unname(draws),
                                         verbose = FALSE)
    expect_identical(colnames(from.response$mean), ids)
})

test_that("conflicting or ambiguous supplied identities are rejected", {
    f <- posterior.path.fixture()
    ids <- paste0("sample-", seq_along(f$field))
    draws <- cbind(f$field, -f$field)
    rownames(draws) <- ids
    reversed <- draws[5:1, , drop = FALSE]
    expect_error(
        lcor.with.posterior(f$adjacency, f$lengths, setNames(f$field, ids),
                            list(marker = reversed), verbose = FALSE),
        "marker.*vertex.*order"
    )
    expect_error(
        lcor.with.posterior(f$adjacency, f$lengths, f$field,
                            list(marker = draws, reversed = reversed), verbose = FALSE),
        "reversed.*vertex.*order"
    )
    for (bad.names in list(c("same", "same"), c("marker", ""), c("marker", NA_character_))) {
        expect_error(
            lcor.with.posterior(f$adjacency, f$lengths, f$field,
                                setNames(list(draws, draws), bad.names), verbose = FALSE),
            "names\\(Z.hat.samples\\)"
        )
    }
    rownames(draws)[1] <- rownames(draws)[2]
    expect_error(
        lcor.with.posterior(f$adjacency, f$lengths, f$field, draws, verbose = FALSE),
        "rownames.*unique"
    )
})

permutation.fixture <- function() {
    a <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
    list(a = a, w = lapply(a, function(v) rep(1, length(v))),
         y = c(1, 4, 2, 3), z = cbind(A = c(3, 1, 2, 5), B = c(2, 4, 1, 3)))
}

test_that("seeded tests restore present and absent RNG states, including errors", {
    f <- permutation.fixture()
    set.seed(83); old <- .Random.seed
    x <- permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = 9)
    expect_identical(.Random.seed, old)
    expect_identical(x$p.value, permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = 9)$p.value)
    expect_error(permutation.test.lcor(f$a, f$w, f$y, f$z, hop.radius = 0), "hop.radius")
    expect_identical(.Random.seed, old)
    permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = 2, seed = NULL)
    expect_false(identical(.Random.seed, old))
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    rm(".Random.seed", envir = .GlobalEnv)
    permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = 2)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("Monte Carlo values agree with an independent edge-cosine calculation", {
    f <- permutation.fixture()
    cosine <- function(z) vapply(seq_along(f$a), function(v) {
        dy <- f$y[f$a[[v]]] - f$y[v]; dz <- z[f$a[[v]]] - z[v]
        den <- sqrt(sum(dy^2) * sum(dz^2))
        if (den <= 1e-10) 0 else sum(dy * dz) / den
    }, numeric(1))
    observed <- apply(f$z, 2, function(z) mean(abs(cosine(z))))
    set.seed(23)
    simulated <- t(replicate(19, apply(f$z[sample.int(4), , drop = FALSE], 2,
                                     function(z) mean(abs(cosine(z))))))
    expected <- (1 + colSums(sweep(simulated, 2, observed, `>=`))) / 20
    x <- permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = 19, seed = 23,
                               return.perm.stats = TRUE)
    expect_equal(x$stat.perm, simulated, tolerance = 1e-12)
    expect_equal(x$p.value, expected)
    expect_equal(x$q.value, p.adjust(expected, "BH"))
    streamed <- permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = 19, seed = 23)
    expect_identical(streamed$p.value, x$p.value)
    expect_null(streamed$stat.perm)
})

test_that("statistic controls and grouped permutations have explicit contracts", {
    f <- permutation.fixture()
    x <- permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = 7, hop.radius = 2,
                               y.diff.type = "logratio", epsilon = 0.25,
                               winsorize.quantile = 0.1, strata = seq_len(4))
    expected <- lcor(f$a, f$w, f$y, f$z, hop.radius = 2, y.diff.type = "logratio",
                     epsilon = 0.25, winsorize.quantile = 0.1)
    expect_equal(unname(x$stat.obs), unname(colMeans(abs(expected))))
    expect_equal(unname(x$p.value), c(1, 1))
    expect_output(print(x, n = 1), "within supplied groups")
    expect_output(print(x), "minimum p-value")
    expect_identical(print(x, n = 0), x)
    expect_error(permutation.test.lcor(f$a, f$w, f$y, f$z, strata = c(1, NA, 2, 2)), "strata")
    expect_error(permutation.test.lcor(f$a, f$w, f$y, f$z, seed = NA), "seed")
    expect_error(permutation.test.lcor(f$a, f$w, f$y, f$z, n.perm = NA), "n.perm")
})

test_that("restricted permutations never cross groups for either permuted field", {
    f <- permutation.fixture(); groups <- list(1:2,3:4)
    for (side in c("y","z")) {
        set.seed(72)
        expected <- t(replicate(9, {
            order <- 1:4
            for (g in groups) order[g] <- g[sample.int(length(g))]
            yy <- if(side == "y") f$y[order] else f$y
            zz <- if(side == "z") f$z[order,,drop=FALSE] else f$z
            colMeans(abs(lcor(f$a,f$w,yy,zz)))
        }))
        result <- permutation.test.lcor(f$a,f$w,f$y,f$z,n.perm=9,seed=72,
                                        permute=side,strata=c(1,1,2,2),return.perm.stats=TRUE)
        expect_equal(result$stat.perm,expected)
    }
})

test_that("permutation print rejects malformed result tables", {
    f <- permutation.fixture()
    x <- permutation.test.lcor(f$a,f$w,f$y,f$z,n.perm=1)
    x$table <- NULL
    expect_error(print(x), "table")
})

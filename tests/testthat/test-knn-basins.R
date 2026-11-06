# Independent dense-distance oracle for small geometry fixtures.
knn.basins.reference.graph <- function(X, k, metric) {
    d <- as.matrix(stats::dist(X, method = if (metric == "L2") "euclidean" else "maximum"))
    diag(d) <- Inf
    index <- matrix(unlist(lapply(seq_len(nrow(X)), function(i) {
        as.integer(order(d[i, ], seq_len(nrow(X)))[seq_len(k)])
    })), nrow = nrow(X), byrow = TRUE)
    directed <- matrix(FALSE, nrow(X), nrow(X))
    directed[cbind(rep(seq_len(nrow(X)), each = k), as.vector(t(index)))] <- TRUE
    list(index = index, adjacency = directed | t(directed), distance = d)
}

# Threshold-connected components give an independent oracle for birth/death
# pairs, without invoking the package's union-find or plateau implementation.
knn.basins.reference.persistence <- function(adjacency, y, peak) {
    reachable <- function(level) {
        seen <- peak
        repeat {
            candidates <- which(colSums(adjacency[seen, , drop = FALSE]) > 0 & y >= level)
            expanded <- union(seen, candidates)
            if (length(expanded) == length(seen)) return(seen)
            seen <- expanded
        }
    }
    for (level in sort(unique(y[y <= y[peak]]), decreasing = TRUE)) {
        region <- reachable(level)
        if (any(y[region] > y[peak])) return(y[peak] - level)
    }
    Inf
}

test_that("exact L2 and Linfty neighbors agree with independent distance matrices", {
    X <- rbind(c(0, 0), c(2, 0), c(1.5, 1.5), c(-1, 2), c(3, 4), c(-3, -2))
    y <- c(0, 2, 1, -1, 3, 0)
    for (metric in c("L2", "Linfty")) {
        for (k in c(1L, 3L, nrow(X) - 1L)) {
            fit <- knn.basins(X, y, k, metric)
            reference <- knn.basins.reference.graph(X, k, metric)
            g <- fit$neighborhood
            expect_identical(g$knn.index, reference$index)
            expect_equal(g$radius, g$knn.distance[, k])
            for (i in seq_len(nrow(X))) {
                expect_identical(g$adj.list[[i]], which(reference$adjacency[i, ]))
                expect_equal(g$edge.length.list[[i]], unname(reference$distance[i, g$adj.list[[i]]]))
                expect_equal(g$knn.distance[i, ], unname(reference$distance[i, g$knn.index[i, ]]))
                expect_false(i %in% g$knn.index[i, ])
            }
        }
    }
    expect_identical(knn.basins(X, y, 1, "L2")$neighborhood$knn.index[1, 1], 2L)
    expect_identical(knn.basins(X, y, 1, "Linfty")$neighborhood$knn.index[1, 1], 3L)
})

test_that("union symmetrization defines extrema as well as persistence", {
    # Row 2 selects row 1, but row 3 selects row 2 and has a higher field.
    fit <- knn.basins(matrix(1:4), c(0, 3, 4, 0), k = 1)
    expect_identical(fit$neighborhood$knn.index[, 1], c(2L, 1L, 2L, 3L))
    expect_identical(fit$max$extrema$row, 3L)
    expect_identical(fit$max$labels, rep(3L, 4))
})

test_that("known peak merges and separate size filtering preserve all stages", {
    X <- matrix(1:9)
    y <- c(0, 4, 3, 1, 2, 3, 2, 1, 0)
    original <- knn.basins(X, y, k = 2)
    tab <- original$max$extrema
    expect_identical(tab$row, c(2L, 6L))
    expect_equal(tab$value, c(4, 3))
    expect_equal(tab$merge.level, c(NA, 1))
    expect_equal(tab$persistence, c(Inf, 2))
    expect_identical(tab$parent.row, c(NA_integer_, 2L))
    expect_identical(tab$initial.size, c(4L, 5L))
    expect_identical(original$max$initial.labels, c(rep(2L, 4), rep(6L, 5)))
    expect_equal(knn.basins(X, y, 2, min.persistence = 2)$max$labels, original$max$labels)

    merged <- knn.basins(X, y, 2, min.persistence = 2.1, min.basin.size = 8)
    expect_identical(merged$max$labels, rep(2L, 9))
    expect_identical(merged$max$initial.labels, original$max$initial.labels)
    expect_identical(merged$max$extrema$merged.size, c(9L, 0L))
    expect_identical(merged$max$extrema$status, c("retained", "merged"))
    expect_identical(merged$max$basins, list(`2` = 1:9))
    filtered <- knn.basins(X, y, 2, min.basin.size = 5)
    expect_identical(filtered$max$labels, c(rep(NA_integer_, 4), rep(6L, 5)))
    expect_identical(filtered$max$persistence.labels, original$max$labels)
    expect_identical(filtered$max$extrema$status, c("filtered_size", "retained"))
    expect_identical(filtered$max$basins, list(`6` = 5:9))
    expect_true(all(is.na(knn.basins(X, y, 2, min.basin.size = 10)$max$labels)))
    expect_length(knn.basins(X, y, 2, min.basin.size = 10)$max$basins, 0)
})

test_that("persistence merging follows a chain of surviving ancestors", {
    X <- matrix(1:13)
    y <- c(0, 10, 6, 2, 5, 9, 8, 7, 7.5, 8, 4, 1, 0)
    fit <- knn.basins(X, y, 2, direction = "max")
    expect_identical(fit$max$extrema$row, c(2L, 6L, 10L))
    expect_equal(fit$max$extrema$persistence, c(Inf, 7, 1))
    expect_identical(fit$max$extrema$parent.row, c(NA_integer_, 2L, 6L))
    partial <- knn.basins(X, y, 2, direction = "max", min.persistence = 2)
    expect_identical(partial$max$extrema$persistence.root, c(2L, 6L, 6L))
    full <- knn.basins(X, y, 2, direction = "max", min.persistence = 8)
    expect_identical(full$max$extrema$persistence.root, rep(2L, 3))
    expect_identical(full$max$labels, rep(2L, 13))
})

test_that("a qualifying extremum can have a singleton basin and be size-filtered", {
    X <- matrix(c(-1.2, -1.1, -1, 0, 1, 1.1, 1.2))
    y <- c(1, 2, 0, 1, 0, 2, 1)
    original <- knn.basins(X, y, 2, direction = "max")
    row <- match(4L, original$max$extrema$row)
    expect_false(is.na(row))
    expect_identical(original$max$extrema$initial.size[row], 1L)
    filtered <- knn.basins(X, y, 2, direction = "max", min.basin.size = 2)
    expect_true(is.na(filtered$max$labels[4]))
    expect_identical(filtered$max$extrema$status[row], "filtered_size")
    expect_true(all(!is.na(filtered$max$labels[-4])))
})

test_that("merging uses the topological elder even across a retained flow basin", {
    # Peak 5 already joined peak 2 at height 10 when peak 8 dies at height 7.
    # Peak 5 remains significant (persistence 9), whereas peak 8 has lifetime 1.
    # Its topological elder is therefore 2, not its adjacent flow basin 5.
    fit <- knn.basins(matrix(1:9), c(0, 20, 12, 10, 19, 9, 7, 8, 0),
                      k = 2, direction = "max", min.persistence = 2)
    expect_identical(fit$max$extrema$row, c(2L, 5L, 8L))
    expect_equal(fit$max$extrema$persistence, c(Inf, 9, 1))
    expect_identical(fit$max$extrema$parent.row, c(NA_integer_, 2L, 2L))
    expect_identical(fit$max$labels, c(rep(2L, 4), rep(5L, 3), rep(2L, 2)))
})

test_that("persistence matches independent threshold connectivity in both metrics", {
    set.seed(8472)
    X <- matrix(rnorm(120), ncol = 4)
    y <- sample(seq_len(nrow(X))) / 10
    for (metric in c("L2", "Linfty")) {
        fit <- knn.basins(X, y, k = 3, metric = metric)
        reference <- knn.basins.reference.graph(X, 3, metric)
        for (direction in c("max", "min")) {
            field <- if (direction == "max") y else -y
            peaks <- which(vapply(seq_along(y), function(i) {
                all(field[i] > field[reference$adjacency[i, ]])
            }, logical(1)))
            observed <- fit[[direction]]$extrema
            expect_setequal(observed$row, peaks)
            expected <- vapply(observed$row, function(i) {
                knn.basins.reference.persistence(reference$adjacency, field, i)
            }, numeric(1))
            expect_equal(observed$persistence, expected)
            expect_equal(sum(observed$initial.size), nrow(X))
            expect_equal(sum(observed$merged.size), nrow(X))
        }
    }
})

test_that("minima are exactly maxima of the negated outcome including filtering", {
    X <- cbind(1:13, c(3, 1, 6, 2, 4, 9, 10, 8, 5, 7, 13, 11, 12))
    y <- c(0, 10, 6, 2, 5, 9, 8, 7, 7, 8, 4, 1, 0)
    for (metric in c("L2", "Linfty")) {
        fit <- knn.basins(X, y, 2, metric, min.persistence = 3, min.basin.size = 3)
        reversed <- knn.basins(X, -y, 2, metric, "max", 3, 3)
        expected <- reversed$max
        expected$extrema$value <- -expected$extrema$value
        expected$extrema$merge.level <- -expected$extrema$merge.level
        expect_equal(fit$min, expected)
        expect_null(reversed$min)
        expect_null(knn.basins(X, y, 2, direction = "min")$max)
    }
})

test_that("plateaus, equal-height peaks, and distance ties are deterministic", {
    X <- matrix(1:9)
    fit <- knn.basins(X, c(0, 3, 3, 1, 2, 4, 4, 0, 0), 2)
    expect_identical(fit$max$extrema$row, c(6L, 2L))
    expect_identical(fit$max$extrema$plateau.rows, I(list(6:7, 2:3)))
    expect_equal(fit$max$extrema$persistence, c(Inf, 2))
    expect_equal(fit$max$extrema$local.gap, c(2, 2))
    expect_identical(fit$max$initial.labels[2], fit$max$initial.labels[3])
    flat <- knn.basins(X, rep(1, 9), 2, min.persistence = Inf)
    for (direction in c("max", "min")) {
        expect_identical(flat[[direction]]$labels, rep(1L, 9))
        expect_identical(flat[[direction]]$extrema$plateau.size, 9L)
        expect_equal(flat[[direction]]$extrema$persistence, Inf)
        expect_true(is.na(flat[[direction]]$extrema$local.gap))
    }
    ties <- knn.basins(matrix(1:7), c(0, 3, 1, 0, 1, 3, 0), 2)
    expect_identical(ties$max$extrema$row, c(2L, 6L))
    expect_identical(ties$max$extrema$parent.row, c(NA_integer_, 2L))
    expect_equal(ties$max$extrema$persistence, c(Inf, 3))
    set.seed(123)
    before <- .Random.seed
    expect_identical(knn.basins(X, rep(1, 9), 2, min.persistence = Inf), flat)
    expect_identical(.Random.seed, before)
})

test_that("disconnected components retain separate essential extrema", {
    X <- matrix(c(0, 1, 2, 100, 101, 102))
    fit <- knn.basins(X, c(0, 2, 1, 5, 4, 3), 2, min.persistence = Inf)
    expect_identical(fit$neighborhood$component, rep(1:2, each = 3))
    expect_identical(fit$neighborhood$n.components, 2L)
    for (direction in c("max", "min")) {
        expect_true(all(fit[[direction]]$extrema$essential))
        expect_true(all(is.infinite(fit[[direction]]$extrema$persistence)))
        expect_equal(length(unique(fit[[direction]]$labels)), 2L)
    }
    filtered <- knn.basins(X, c(0, 2, 1, 5, 4, 3), 2,
                           min.persistence = Inf, min.basin.size = 4)
    expect_true(all(is.na(filtered$max$labels)))
    expect_true(all(is.na(filtered$min$labels)))
})

test_that("duplicate rows and zero-length edges are supported without perturbation", {
    X <- matrix(c(0, 0, 1, 2, 2))
    rownames(X) <- letters[1:5]
    fit <- knn.basins(X, c(2, 1, 3, 0, 0), 2)
    expect_identical(fit$row.names, letters[1:5])
    expect_equal(fit$neighborhood$knn.distance[1, 1], 0)
    expect_identical(fit$max$labels, rep(3L, 5))
    expect_identical(fit$min$labels[4], fit$min$labels[5])
    expect_true(all(is.finite(fit$neighborhood$radius)))
    expect_identical(fit$parameters$projection, "none")
    expect_identical(fit$parameters$scaling, "none")
})

test_that("inputs and numeric distance range are validated", {
    X <- matrix(1:6, ncol = 2)
    y <- 1:3
    for (bad in list(as.data.frame(X), matrix(1, 1, 2), matrix(numeric(), 3, 0),
                     matrix(c(NA, 1, 2)), matrix(c(Inf, 1, 2)), X + 1i)) {
        expect_error(knn.basins(bad, y, 2), "X")
    }
    for (bad in list(1:2, c(1, NA, 3), c(1, Inf, 3), matrix(y), letters[1:3], y + 1i)) {
        expect_error(knn.basins(X, bad, 2), "y")
    }
    for (bad in list(0, 3, 1.5, NA, Inf, c(1, 2), TRUE)) {
        expect_error(knn.basins(X, y, bad), "k")
    }
    for (bad in list(-1, NA, NaN, c(0, 1), "1")) {
        expect_error(knn.basins(X, y, 2, min.persistence = bad), "min.persistence")
    }
    for (bad in list(0, -1, 1.5, NA, Inf, c(1, 2))) {
        expect_error(knn.basins(X, y, 2, min.basin.size = bad), "min.basin.size")
    }
    for (bad in list(NA, 0, 1, "FALSE", c(TRUE, FALSE), logical(), NULL)) {
        expect_error(knn.basins(X, y, 2, symmetrize = bad), "symmetrize")
    }
    expect_error(knn.basins(X, y, 2, metric = "linf.simplex"), "metric")
    expect_error(knn.basins(X, y, 2, direction = "up"), "direction")
    expect_error(knn.basins(matrix(c(-1e308, 1e308)), c(0, 1), 1), "rescale X")
    expect_error(knn.basins(X, c(-1e308, 0, 1e308), 2), "rescale y")
    huge <- knn.basins(matrix(c(0, 1e200)), c(0, 1), 1)
    tiny <- knn.basins(matrix(c(0, 1e-200)), c(0, 1), 1)
    expect_equal(huge$neighborhood$radius / 1e200, c(1, 1))
    expect_equal(tiny$neighborhood$radius / 1e-200, c(1, 1))
    expect_output(print(huge), "max: 1 original extrema")
    expect_error(print(structure(list(), class = "gflow_knn_basins")), "missing fields")
    broken <- huge
    broken$max$labels <- integer()
    expect_error(print(broken), "malformed max")
})

# Independent directed oracle: enumerate field thresholds and traverse only
# outgoing edges. Neither widest-path code nor the package's flow is called.
knn.basins.reference.escape <- function(adjacency, field, peaks, peak) {
    for (level in sort(unique(field[field <= field[peak]]), decreasing = TRUE)) {
        seen <- peak
        repeat {
            add <- which(colSums(adjacency[seen, , drop = FALSE]) > 0 & field >= level)
            expanded <- union(seen, add)
            if (length(expanded) == length(seen)) break
            seen <- expanded
        }
        targets <- intersect(seen, peaks[field[peaks] > field[peak]])
        if (length(targets)) {
            parent <- targets[order(-field[targets], targets)][1]
            return(list(level = level, parent = parent, persistence = field[peak] - level))
        }
    }
    list(level = NA_real_, parent = NA_integer_, persistence = Inf)
}

test_that("asymmetric neighborhoods test exactly k outgoing choices", {
    X <- matrix(c(0, 1, 3.5))
    y <- c(0, 5, 10)
    directed <- knn.basins(X, y, 1, symmetrize = FALSE)
    symmetric <- knn.basins(X, y, 1)
    expect_identical(symmetric$max$extrema$row, 3L)
    expect_identical(directed$max$extrema$row, c(3L, 2L))
    expect_identical(directed$max$labels, c(2L, 2L, 3L))
    expect_identical(directed$max$next.vertex, c(2L, NA_integer_, NA_integer_))
    expect_identical(directed$min$labels, rep(1L, 3))
    expect_equal(lengths(directed$neighborhood$adj.list), rep(1L, 3))
    expect_identical(directed$neighborhood$adj.list[[2]], 1L)
    expect_equal(directed$max$extrema$local.gap, c(5, 5))
    # The lower peak cannot follow the incoming edge from the higher peak.
    expect_equal(directed$neighborhood$n.components, 1L)
    expect_true(all(directed$max$extrema$essential))
    expect_true(all(is.infinite(directed$max$extrema$persistence)))
    expect_identical(directed$parameters$symmetrization, "none")
    expect_identical(directed$parameters$persistence, "directed_escape")
    expect_output(print(directed), "directed escape prominence")
    explicit <- knn.basins(X, y, 1, symmetrize = TRUE)
    expect_identical(explicit$max, symmetric$max)
    expect_identical(explicit$min, symmetric$min)
    expect_identical(explicit$neighborhood, symmetric$neighborhood)
})

test_that("directed flow follows the graph CLOSEST shortest improving edge rule", {
    # Row 1 can ascend one unit to the modest peak or two units to a much
    # higher peak. The short step determines its initial basin.
    X <- matrix(c(0, 1, 2.2, -2, -2.1))
    y <- c(0, 2, 1, 100, 99)
    fit <- knn.basins(X, y, 2, symmetrize = FALSE)
    expect_identical(fit$max$next.vertex[1], 2L)
    expect_identical(fit$max$initial.labels, c(2L, 2L, 2L, 4L, 4L))
    # Direct parity with the selector used by canonical plateau trajectory
    # flow, applied to singleton plateaus and no edge-length cutoff.
    graph <- fit$neighborhood
    plateaus <- list(plateau.of.vertex = seq_along(y), level = y)
    edges <- gflow:::.trajectory.plateau.boundary.edges(
        plateaus, graph$adj.list, graph$edge.length.list
    )
    for (direction in c("max", "min")) for (i in seq_along(y)) {
        selected <- gflow:::.trajectory.plateau.outgoing.edge(
            i, direction, plateaus, edges, Inf, "allow"
        )$edge
        expected <- if (is.null(selected)) NA_integer_ else selected$target.vertex
        expect_identical(fit[[direction]]$next.vertex[i], expected)
    }
    ties <- knn.basins(matrix(c(-1, 0, 1)), c(1, 0, 100), 2, symmetrize = FALSE)
    expect_identical(ties$max$next.vertex[2], 1L)
})

test_that("directed extrema and escape levels match independent threshold reachability", {
    set.seed(19824)
    X <- matrix(rnorm(54), ncol = 3)
    # Repeated field values exercise weak extrema and equal escape targets.
    y <- sample(rep(1:9, each = 2))
    for (metric in c("L2", "Linfty")) for (k in c(1L, 3L, 8L)) {
        fit <- knn.basins(X, y, k, metric, symmetrize = FALSE)
        reference <- knn.basins.reference.graph(X, k, metric)
        adjacency <- matrix(FALSE, nrow(X), nrow(X))
        adjacency[cbind(rep(seq_len(nrow(X)), each = k),
                         as.vector(t(reference$index)))] <- TRUE
        expect_identical(fit$neighborhood$knn.index, reference$index)
        for (i in seq_len(nrow(X))) {
            expect_identical(fit$neighborhood$adj.list[[i]], which(adjacency[i, ]))
            expect_equal(fit$neighborhood$edge.length.list[[i]],
                         unname(reference$distance[i, adjacency[i, ]]))
        }
        for (direction in c("max", "min")) {
            field <- if (direction == "max") y else -y
            peaks <- which(vapply(seq_along(y), function(i)
                all(field[i] >= field[adjacency[i, ]]), logical(1)))
            tab <- fit[[direction]]$extrema
            expect_setequal(tab$row, peaks)
            for (i in seq_len(nrow(tab))) {
                oracle <- knn.basins.reference.escape(adjacency, field, peaks, tab$row[i])
                expect_equal(tab$persistence[i], oracle$persistence)
                expect_equal(tab$merge.level[i], if (direction == "max") oracle$level else -oracle$level)
                expect_identical(tab$parent.row[i], oracle$parent)
            }
            next.vertex <- fit[[direction]]$next.vertex
            for (i in seq_along(y)) {
                better <- which(adjacency[i, ] & field > field[i])
                expected <- if (!length(better)) NA_integer_ else
                    better[order(reference$distance[i, better], better)][1]
                expect_identical(next.vertex[i], expected)
                terminal <- i
                while (!is.na(next.vertex[terminal])) terminal <- next.vertex[terminal]
                expect_identical(fit[[direction]]$initial.labels[i], terminal)
            }
            expect_equal(sum(tab$initial.size), nrow(X))
            expect_false(anyNA(fit[[direction]]$initial.labels))
        }
    }
})

test_that("directed ties keep the per-row test and do not pool neighborhoods", {
    # Row 1 sees an equal row and stops. Row 2 sees row 3 and ascends.
    fit <- knn.basins(matrix(c(0, 2, 3)), c(5, 5, 10), 1, symmetrize = FALSE)
    expect_identical(fit$max$extrema$row, c(3L, 1L))
    expect_equal(fit$max$extrema$persistence, c(Inf, 0))
    expect_equal(fit$max$extrema$merge.level, c(NA, 5))
    expect_identical(fit$max$extrema$parent.row, c(NA_integer_, 3L))
    expect_identical(fit$max$initial.labels, c(1L, 3L, 3L))
    expect_identical(fit$max$extrema$plateau.size, c(1L, 1L))
    merged <- knn.basins(matrix(c(0, 2, 3)), c(5, 5, 10), 1,
                         symmetrize = FALSE, min.persistence = 0.1)
    expect_identical(merged$max$labels, rep(3L, 3))
    expect_identical(merged$max$next.vertex, fit$max$next.vertex)
    flat <- knn.basins(matrix(1:6), rep(7, 6), 2,
                       symmetrize = FALSE, min.persistence = Inf)
    for (d in c("max", "min")) {
        expect_identical(flat[[d]]$initial.labels, 1:6)
        expect_true(all(flat[[d]]$extrema$essential))
        expect_equal(flat[[d]]$extrema$local.gap, rep(0, 6))
        expect_true(all(is.na(flat[[d]]$next.vertex)))
    }
    # Equal-height peaks are not strictly better escape targets.
    equal <- knn.basins(matrix(1:5), c(0, 3, 1, 3, 0), 2, symmetrize = FALSE)
    expect_identical(equal$max$extrema$row, c(2L, 4L))
    expect_true(all(equal$max$extrema$essential))
})

test_that("directed prominence merging precedes independent size filtering", {
    X <- matrix(1:9); y <- c(0, 4, 3, 1, 2, 3, 2, 1, 0)
    fit <- knn.basins(X, y, 2, symmetrize = FALSE, min.persistence = 2)
    expect_identical(fit$max$extrema$row, c(2L, 6L))
    expect_equal(fit$max$extrema$persistence, c(Inf, 2))
    expect_identical(fit$max$labels, c(rep(2L, 4), rep(6L, 5)))
    merged <- knn.basins(X, y, 2, symmetrize = FALSE,
                         min.persistence = 2.1, min.basin.size = 8)
    expect_identical(merged$max$labels, rep(2L, 9))
    expect_identical(merged$max$initial.labels, fit$max$initial.labels)
    expect_identical(merged$max$next.vertex, fit$max$next.vertex)
    filtered <- knn.basins(X, y, 2, symmetrize = FALSE, min.basin.size = 5)
    expect_identical(filtered$max$labels, c(rep(NA_integer_, 4), rep(6L, 5)))
    expect_identical(filtered$max$persistence.labels, fit$max$labels)
    chain <- knn.basins(matrix(1:13), c(0,10,6,2,5,9,8,7,7.5,8,4,1,0), 2,
                        symmetrize = FALSE, min.persistence = 8)
    expect_identical(chain$max$extrema$parent.row, c(NA_integer_, 2L, 6L))
    expect_identical(chain$max$labels, rep(2L, 13))
})

test_that("directed minima reverse maxima including duplicate coordinates and ties", {
    X <- cbind(c(0,0,1,2,2,4,7), c(0,0,3,1,1,0,4))
    y <- c(0,2,2,1,3,0,3)
    for (metric in c("L2", "Linfty")) {
        fit <- knn.basins(X, y, 2, metric, symmetrize = FALSE,
                          min.persistence = 1.5, min.basin.size = 2)
        inverse <- knn.basins(X, -y, 2, metric, direction = "max", symmetrize = FALSE,
                              min.persistence = 1.5, min.basin.size = 2)$max
        inverse$extrema$value <- -inverse$extrema$value
        inverse$extrema$merge.level <- -inverse$extrema$merge.level
        expect_identical(fit$min, inverse)
        expect_true(all(is.finite(fit$neighborhood$knn.distance)))
    }
})

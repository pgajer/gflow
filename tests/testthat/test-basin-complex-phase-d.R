phase.d.graph <- function(n, edges) {
    adj.list <- vector("list", n)
    edge.length.list <- vector("list", n)
    for (i in seq_len(nrow(edges))) {
        from <- as.integer(edges[i, 1L])
        to <- as.integer(edges[i, 2L])
        edge.length <- edges[i, 3L]
        adj.list[[from]] <- c(adj.list[[from]], to)
        adj.list[[to]] <- c(adj.list[[to]], from)
        edge.length.list[[from]] <- c(
            edge.length.list[[from]],
            edge.length
        )
        edge.length.list[[to]] <- c(edge.length.list[[to]], edge.length)
    }
    for (vertex in seq_len(n)) {
        if (is.null(adj.list[[vertex]])) {
            adj.list[[vertex]] <- integer()
            edge.length.list[[vertex]] <- numeric()
        }
    }
    list(
        adj.list = adj.list,
        edge.length.list = edge.length.list
    )
}

phase.d.grid.fixture <- function() {
    edges <- matrix(
        c(
            1, 2, 1, 2, 3, 1,
            4, 5, 1, 5, 6, 1,
            7, 8, 1, 8, 9, 1,
            1, 4, 1, 4, 7, 1,
            2, 5, 1, 5, 8, 1,
            3, 6, 1, 6, 9, 1
        ),
        ncol = 3L,
        byrow = TRUE
    )
    c(
        phase.d.graph(9L, edges),
        list(
            field = c(4, 4, 1, 2, 0, 1, 1, 2, 3),
            mass = c(2, 1, 1, 2, 3, 1, 1, 2, 2)
        )
    )
}

test_that("merge tree contracts connected exact plateaus on a 2D grid", {
    fixture <- phase.d.grid.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "superlevel_merge_tree",
        direction = "max",
        vertex.mass = fixture$mass
    )

    expect_identical(object$status, "ok")
    expect_silent(gflow:::.validate.basin.complex(object))
    plateau.row <- which(object$extrema$representative.vertex == 1L)
    expect_length(plateau.row, 1L)
    expect_identical(
        object$extrema$plateau.vertices[[plateau.row]],
        c(1L, 2L)
    )
    expect_identical(
        object$extrema$n.plateau.vertices[[plateau.row]],
        2L
    )
    expect_true(all(object$basin.table$persistence >= 0))
    expect_equal(nrow(object$assignment), 9L)
    expect_true(all(object$assignment$assignment.status == "assigned"))
})

test_that("maximum and minimum trees use opposite persistence signs", {
    fixture <- phase.d.grid.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "superlevel_merge_tree",
        direction = "both"
    )

    expect_true(all(object$basin.table$persistence >= 0))
    max.rows <- object$basin.table$type == "max"
    min.rows <- object$basin.table$type == "min"
    expect_equal(
        object$basin.table$persistence[max.rows],
        object$basin.table$birth.level[max.rows] -
            object$basin.table$death.level[max.rows]
    )
    expect_equal(
        object$basin.table$persistence[min.rows],
        object$basin.table$death.level[min.rows] -
            object$basin.table$birth.level[min.rows]
    )
    expect_equal(
        as.integer(table(object$assignment$direction)),
        c(9L, 9L)
    )
})

test_that("elder ties resolve by representative vertex", {
    graph <- phase.d.graph(
        3L,
        matrix(c(1, 2, 1, 2, 3, 1), ncol = 3L, byrow = TRUE)
    )
    object <- create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(3, 1, 3),
        method = "superlevel_merge_tree",
        direction = "max"
    )

    expect_equal(nrow(object$merge.table), 1L)
    expect_identical(object$merge.table$event.status, "elder_tie")
    expect_identical(
        object$merge.table$surviving.basin.id,
        "basin_max_v00000001"
    )
    expect_identical(
        object$merge.table$losing.basin.id,
        "basin_max_v00000003"
    )
    child <- object$basin.table$basin.id == "basin_max_v00000003"
    expect_identical(
        object$basin.table$parent.basin.id[child],
        "basin_max_v00000001"
    )
})

test_that("one saddle can record multiple deterministic merges", {
    graph <- phase.d.graph(
        4L,
        matrix(
            c(1, 4, 1, 2, 4, 1, 3, 4, 1),
            ncol = 3L,
            byrow = TRUE
        )
    )
    object <- create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(4, 3, 2, 0),
        method = "superlevel_merge_tree",
        direction = "max"
    )

    expect_equal(nrow(object$merge.table), 2L)
    expect_true(all(
        object$merge.table$merge.plateau.id == "plateau_00000004"
    ))
    expect_true(all(
        object$merge.table$surviving.basin.id ==
            "basin_max_v00000001"
    ))
})

test_that("disconnected components and isolates have independent roots", {
    graph <- phase.d.graph(
        5L,
        matrix(c(1, 2, 1, 3, 4, 1), ncol = 3L, byrow = TRUE)
    )
    object <- create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(3, 0, 3, 0, 2),
        method = "superlevel_merge_tree",
        direction = "max"
    )

    expect_equal(nrow(object$basin.table), 3L)
    expect_true(all(is.na(object$basin.table$parent.basin.id)))
    expect_equal(
        sort(object$basin.table$persistence),
        c(0, 3, 3)
    )
    expect_identical(
        object$assignment$basin.id[object$assignment$vertex == 5L],
        "basin_max_v00000005"
    )
})

test_that("tree output is invariant to adjacency row order", {
    fixture <- phase.d.grid.fixture()
    reversed.adj <- lapply(fixture$adj.list, rev)
    reversed.length <- lapply(fixture$edge.length.list, rev)
    first <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "superlevel_merge_tree"
    )
    second <- create.basin.complex(
        reversed.adj,
        reversed.length,
        fixture$field,
        method = "superlevel_merge_tree"
    )

    expect_identical(first$extrema, second$extrema)
    expect_identical(first$basin.table, second$basin.table)
    expect_identical(first$membership, second$membership)
    expect_identical(first$assignment, second$assignment)
    expect_identical(first$merge.table, second$merge.table)
})

test_that("hierarchical support and primary mass retain distinct meanings", {
    fixture <- phase.d.grid.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "superlevel_merge_tree",
        direction = "both",
        vertex.mass = fixture$mass
    )

    weight.sum <- aggregate(
        membership.weight ~ vertex + direction,
        object$membership,
        sum
    )
    expect_equal(weight.sum$membership.weight, rep(1, nrow(weight.sum)))
    primary.by.direction <- tapply(
        object$basin.table$primary.support.mass,
        object$basin.table$type,
        sum
    )
    expect_equal(as.numeric(primary.by.direction), c(1, 1))
    expect_true(any(
        object$basin.table$raw.support.size >
            object$basin.table$primary.support.size
    ))
})

test_that("merge-tree validation rejects invalid hierarchy and persistence", {
    fixture <- phase.d.grid.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "superlevel_merge_tree"
    )

    bad.persistence <- object
    bad.persistence$basin.table$persistence[[1L]] <- -1
    expect_error(
        gflow:::.validate.basin.complex(bad.persistence),
        class = "gflow_basin_schema_error"
    )

    if (nrow(object$merge.table) > 0L) {
        bad.parent <- object
        child <- match(
            bad.parent$merge.table$losing.basin.id[[1L]],
            bad.parent$basin.table$basin.id
        )
        bad.parent$basin.table$parent.basin.id[[child]] <- NA_character_
        expect_error(
            gflow:::.validate.basin.complex(bad.parent),
            class = "gflow_basin_schema_error"
        )
    }
})

merge.tree.test.graph <- function(n, edges) {
    adjacency <- edge.lengths <- vector("list", n)
    for (index in seq_len(nrow(edges))) {
        from <- edges[index, 1L]
        to <- edges[index, 2L]
        length <- edges[index, 3L]
        adjacency[[from]] <- c(adjacency[[from]], to)
        adjacency[[to]] <- c(adjacency[[to]], from)
        edge.lengths[[from]] <- c(edge.lengths[[from]], length)
        edge.lengths[[to]] <- c(edge.lengths[[to]], length)
    }
    list(adj.list = adjacency, edge.length.list = edge.lengths)
}

merge.tree.test.fixture <- function(direction = "max") {
    graph <- merge.tree.test.graph(
        4L,
        matrix(
            c(1, 4, 1, 2, 4, 1, 3, 4, 1),
            ncol = 3L,
            byrow = TRUE
        )
    )
    field <- c(5, 4, 3, 1)
    if (direction == "min") {
        field <- -field
    }
    create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        field,
        method = "superlevel_merge_tree",
        direction = direction,
        vertex.mass = c(4, 3, 2, 1),
        vertex.id = paste0("sample_", seq_len(4L))
    )
}

merge.tree.test.relational.corruptions <- function(tree, direction) {
    branch.rows <- which(tree$basin.table$type == direction)
    child.row <- branch.rows[
        !is.na(tree$basin.table$parent.basin.id[branch.rows])
    ][[1L]]
    root.row <- branch.rows[
        is.na(tree$basin.table$parent.basin.id[branch.rows])
    ][[1L]]
    child.id <- tree$basin.table$basin.id[[child.row]]
    parent.id <- tree$basin.table$parent.basin.id[[child.row]]
    event.row <- which(
        tree$merge.table$direction == direction &
            tree$merge.table$losing.basin.id == child.id
    )[[1L]]
    alternate.parent <- setdiff(
        tree$basin.table$basin.id[branch.rows],
        c(child.id, parent.id)
    )[[1L]]

    mutate.tree <- function(branch = NULL, event = NULL) {
        corrupted <- tree
        if (!is.null(branch)) {
            corrupted$basin.table[child.row, names(branch)] <-
                unname(branch)
        }
        if (!is.null(event)) {
            corrupted$merge.table[event.row, names(event)] <-
                unname(event)
        }
        corrupted
    }

    parent.disagreement <- mutate.tree(
        branch = c(parent.basin.id = alternate.parent)
    )
    survivor.disagreement <- mutate.tree(
        event = c(surviving.basin.id = alternate.parent)
    )
    event.merge.disagreement <- mutate.tree(
        event = c(
            merge.level =
                tree$merge.table$merge.level[[event.row]] + 0.25
        )
    )
    event.birth.disagreement <- mutate.tree(
        event = c(
            birth.level =
                tree$merge.table$birth.level[[event.row]] + 0.25
        )
    )
    event.death.disagreement <- mutate.tree(
        event = c(
            death.level =
                tree$merge.table$death.level[[event.row]] + 0.25
        )
    )
    event.persistence.disagreement <- mutate.tree(
        event = c(
            persistence =
                tree$merge.table$persistence[[event.row]] + 0.25
        )
    )

    coherent.death <- tree$basin.table$death.level[[child.row]] + 0.25
    branch.death.relation <- mutate.tree(
        branch = c(death.level = coherent.death),
        event = c(
            death.level = coherent.death,
            merge.level = coherent.death
        )
    )
    coherent.persistence <-
        tree$basin.table$persistence[[child.row]] + 0.25
    branch.persistence.relation <- mutate.tree(
        branch = c(persistence = coherent.persistence),
        event = c(persistence = coherent.persistence)
    )

    root.convention <- tree
    root.death <- tree$basin.table$death.level[[root.row]] + 0.25
    root.convention$basin.table$death.level[[root.row]] <- root.death
    root.convention$basin.table$persistence[[root.row]] <-
        if (direction == "max") {
            tree$basin.table$birth.level[[root.row]] - root.death
        } else {
            root.death - tree$basin.table$birth.level[[root.row]]
        }

    list(
        branch_parent = parent.disagreement,
        event_survivor = survivor.disagreement,
        event_merge = event.merge.disagreement,
        event_birth = event.birth.disagreement,
        event_death = event.death.disagreement,
        event_persistence = event.persistence.disagreement,
        branch_death_relation = branch.death.relation,
        branch_persistence_relation = branch.persistence.relation,
        component_root_convention = root.convention
    )
}

test_that("merge-tree getter returns a complete canonical object", {
    complex <- merge.tree.test.fixture()
    tree <- get.basin.merge.tree(complex)

    expect_s3_class(tree, "basin.merge.tree")
    expect_identical(tree$method, "superlevel_merge_tree")
    expect_identical(tree$n.vertices, 4L)
    expect_identical(tree$graph.input$adj.list, complex$graph.input$adj.list)
    expect_identical(tree$field, complex$field)
    expect_identical(tree$basin.table, complex$basin.table)
    expect_identical(tree$merge.table, complex$merge.table)
    expect_identical(tree$membership, complex$membership)
    expect_identical(tree$assignment, complex$assignment)
    expect_identical(tree$plateau.graph, complex$raw$plateau.graph)
    expect_silent(gflow:::.validate.basin.merge.tree(tree))
    expect_true(all(c(
        "plot.basin.merge.tree",
        "cut.basin.merge.tree",
        "as.dendrogram.basin.merge.tree"
    ) %in% getNamespaceExports("gflow")))
})

test_that("dendrogram conversion preserves exact hierarchy metadata", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())
    dendrogram <- as.dendrogram.basin.merge.tree(
        tree,
        labels = c(
            basin_max_v00000001 = "M1",
            basin_max_v00000002 = "M2",
            basin_max_v00000003 = "M3"
        )
    )

    expect_s3_class(dendrogram, "dendrogram")
    expect_setequal(labels(dendrogram), c("M1", "M2", "M3"))
    events <- attr(dendrogram, "basin.merge.tree.events")
    expect_equal(events$merge.level, c(1, 1))
    expect_setequal(
        events$losing.basin.id,
        c("basin_max_v00000002", "basin_max_v00000003")
    )
    expect_true(all(diff(events$transformed.height) >= 0))
    expect_identical(attr(dendrogram, "direction"), "max")
    expect_identical(attr(dendrogram, "graph.component"), 1L)
})

test_that("equal-height descendant events precede their ancestors", {
    graph <- merge.tree.test.graph(
        5L,
        matrix(
            c(1, 4, 1, 4, 3, 1, 3, 2, 1, 2, 5, 1),
            ncol = 3L,
            byrow = TRUE
        )
    )
    complex <- create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(5, 1, 4, 1, 3),
        method = "superlevel_merge_tree",
        direction = "max"
    )
    tree <- get.basin.merge.tree(complex)
    dendrogram <- as.dendrogram(tree)
    events <- attr(dendrogram, "basin.merge.tree.events")

    expect_identical(
        events$losing.basin.id,
        c("basin_max_v00000005", "basin_max_v00000003")
    )
    expect_equal(events$merge.level, c(1, 1))
    expect_setequal(
        labels(dendrogram),
        tree$basin.table$basin.id
    )
})

test_that("exact cuts include saddle plateaus at equality", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())

    above <- cut.basin.merge.tree(tree, height = 2)
    expect_identical(above$n.active.vertices, 3L)
    expect_identical(above$n.components, 3L)
    expect_equal(sort(above$components$support.mass), c(0.2, 0.3, 0.4))

    at.saddle <- cut.basin.merge.tree(tree, height = 1)
    expect_identical(at.saddle$n.active.vertices, 4L)
    expect_identical(at.saddle$n.components, 1L)
    expect_identical(
        at.saddle$components$basin.id,
        "basin_max_v00000001"
    )
    expect_identical(
        at.saddle$components$vertices[[1L]],
        seq_len(4L)
    )
    expect_identical(
        at.saddle$components$vertex.ids[[1L]],
        paste0("sample_", seq_len(4L))
    )
    expect_equal(at.saddle$components$support.mass, 1)
})

test_that("minimum-tree cuts use sublevel sets", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture("min"))
    cut <- cut.basin.merge.tree(tree, height = -2, direction = "min")

    expect_identical(cut$relation, "<=")
    expect_identical(cut$n.active.vertices, 3L)
    expect_identical(cut$n.components, 3L)
})

test_that("merge-tree plot renders all optional annotations", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())
    plot.file <- tempfile(fileext = ".pdf")
    grDevices::pdf(plot.file, width = 10, height = 8)
    result <- plot.basin.merge.tree(
        tree,
        label = "extremum.vertex",
        show.mass = TRUE,
        show.support = TRUE,
        show.barcode.guides = TRUE,
        show.barcode.birth.labels = TRUE,
        show.barcode.parent.labels = TRUE,
        height = 2
    )
    grDevices::dev.off()

    expect_gt(file.info(plot.file)$size, 0)
    expect_identical(nrow(result$branches), 3L)
    expect_identical(sort(result$layout$order), 1:3)
    expect_identical(length(result$coordinates$leaf.x), 3L)
    parent.index <- match(
        result$layout$events$surviving.basin.id,
        result$branches$basin.id
    )
    expect_equal(
        result$coordinates$node.x,
        result$coordinates$leaf.x[parent.index]
    )
    unlink(plot.file)
})

test_that("standard S3 generics dispatch to merge-tree methods", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())

    expect_s3_class(as.dendrogram(tree), "dendrogram")
    expect_s3_class(cut(tree, height = 2), "basin.merge.tree.cut")

    plot.file <- tempfile(fileext = ".pdf")
    grDevices::pdf(plot.file, width = 8, height = 6)
    expect_invisible(plot(
        tree,
        type = "barcode",
        label = "extremum.vertex"
    ))
    grDevices::dev.off()
    expect_gt(file.info(plot.file)$size, 0)
    unlink(plot.file)
})

test_that("dendrogram and plot require a selected component for forests", {
    graph <- merge.tree.test.graph(
        4L,
        matrix(c(1, 2, 1, 3, 4, 1), ncol = 3L, byrow = TRUE)
    )
    complex <- create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(3, 0, 2, 0),
        method = "superlevel_merge_tree",
        direction = "max"
    )
    tree <- get.basin.merge.tree(complex)

    expect_error(
        as.dendrogram.basin.merge.tree(tree),
        "merge forest"
    )
    expect_s3_class(
        as.dendrogram.basin.merge.tree(tree, component = 1),
        "dendrogram"
    )
    cut <- cut.basin.merge.tree(tree, height = 0)
    expect_identical(cut$n.components, 2L)
})

test_that("public layout accessor returns complete canonical plot inputs", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())
    layout <- get.basin.merge.tree.layout(tree)

    expect_s3_class(layout, "basin.merge.tree.layout")
    expect_identical(layout$validation.status, "ok")
    expect_identical(layout$direction, "max")
    expect_identical(layout$component, 1L)
    expect_identical(layout$requested.basin.ids, layout$basin.ids)
    expect_length(layout$closure.added.ids, 0L)
    expect_identical(
        layout$component.root.basin.id,
        "basin_max_v00000001"
    )
    expect_identical(layout$leaf.order, layout$basin.ids[layout$order])
    expect_setequal(
        layout$events$losing.basin.id,
        setdiff(
            layout$basin.ids,
            layout$component.root.basin.id
        )
    )
    expect_equal(
        layout$branches$birth.level,
        tree$basin.table$birth.level
    )
    expect_equal(
        layout$branches$death.level,
        tree$basin.table$death.level
    )
    expect_equal(
        layout$branches$persistence,
        tree$basin.table$persistence
    )
    expect_identical(
        layout$branches$parent.basin.id,
        tree$basin.table$parent.basin.id
    )
    canonical.events <- tree$merge.table[
        match(layout$events$event.id, tree$merge.table$event.id),
        ,
        drop = FALSE
    ]
    row.names(canonical.events) <- seq_len(nrow(canonical.events))
    expect_identical(layout$events, canonical.events)
    expect_identical(
        layout$coordinates$branches$basin.id,
        layout$branches$basin.id
    )
    expect_identical(
        layout$coordinates$events$event.id,
        layout$events$event.id
    )
    expect_equal(
        layout$coordinates$events$merge.level,
        layout$events$merge.level
    )
    expect_true(
        "get.basin.merge.tree.layout" %in% getNamespaceExports("gflow")
    )
})

test_that("restricted layouts require or construct canonical closure", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())
    root.id <- "basin_max_v00000001"
    child.id <- "basin_max_v00000003"

    expect_error(
        get.basin.merge.tree.layout(tree, basin.ids = child.id),
        "not ancestor-closed"
    )
    closed <- get.basin.merge.tree.layout(
        tree,
        basin.ids = child.id,
        close.ancestors = TRUE
    )
    expect_identical(closed$requested.basin.ids, child.id)
    expect_identical(closed$closure.added.ids, root.id)
    expect_identical(closed$basin.ids, c(root.id, child.id))
    expect_identical(closed$component.root.basin.id, root.id)
    expect_identical(nrow(closed$events), 1L)
    expect_identical(
        closed$events$losing.basin.id,
        child.id
    )
    expect_identical(
        closed$events$surviving.basin.id,
        root.id
    )

    root.only <- get.basin.merge.tree.layout(
        tree,
        basin.ids = root.id
    )
    expect_identical(root.only$basin.ids, root.id)
    expect_identical(nrow(root.only$events), 0L)
    expect_identical(dim(root.only$merge), c(0L, 2L))
    expect_identical(root.only$order, 1L)
    expect_equal(root.only$coordinates$branches$x, 1)
    expect_identical(nrow(root.only$coordinates$events), 0L)
})

test_that("restricted layout preserves complete canonical leaf order", {
    graph <- merge.tree.test.graph(
        6L,
        matrix(
            c(
                1, 4, 1,
                2, 4, 1,
                3, 5, 1,
                4, 5, 1,
                5, 6, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    tree <- get.basin.merge.tree(create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(8, 7, 6, 2, 1, 0),
        method = "superlevel_merge_tree",
        direction = "max"
    ))
    complete <- get.basin.merge.tree.layout(tree)
    child.ids <- tree$basin.table$basin.id[
        !is.na(tree$basin.table$parent.basin.id)
    ]
    restricted <- get.basin.merge.tree.layout(
        tree,
        basin.ids = child.ids[[length(child.ids)]],
        close.ancestors = TRUE
    )

    expect_identical(
        restricted$leaf.order,
        complete$leaf.order[
            complete$leaf.order %in% restricted$basin.ids
        ]
    )
    expect_equal(
        restricted$branches$birth.level,
        tree$basin.table$birth.level[
            match(restricted$basin.ids, tree$basin.table$basin.id)
        ]
    )
    expect_equal(
        restricted$branches$death.level,
        tree$basin.table$death.level[
            match(restricted$basin.ids, tree$basin.table$basin.id)
        ]
    )
    expect_equal(
        restricted$events$merge.level,
        tree$merge.table$merge.level[
            match(
                restricted$events$event.id,
                tree$merge.table$event.id
            )
        ]
    )
})

test_that("layout accessor rejects invalid canonical selections", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())

    expect_error(
        get.basin.merge.tree.layout(tree, basin.ids = character()),
        "nonempty character vector"
    )
    expect_error(
        get.basin.merge.tree.layout(
            tree,
            basin.ids = rep(tree$basin.table$basin.id[[1L]], 2L)
        ),
        "unique"
    )
    expect_error(
        get.basin.merge.tree.layout(tree, basin.ids = "not_a_basin"),
        "Unknown canonical basin id"
    )

    both <- get.basin.merge.tree(create.basin.complex(
        tree$graph.input$adj.list,
        tree$graph.input$edge.length.list,
        tree$field$construction.values,
        method = "superlevel_merge_tree",
        direction = "both"
    ))
    minimum.id <- both$basin.table$basin.id[
        both$basin.table$type == "min"
    ][[1L]]
    expect_error(
        get.basin.merge.tree.layout(
            both,
            direction = "max",
            basin.ids = minimum.id
        ),
        "different direction"
    )

    graph <- merge.tree.test.graph(
        4L,
        matrix(c(1, 2, 1, 3, 4, 1), ncol = 3L, byrow = TRUE)
    )
    forest <- get.basin.merge.tree(create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(3, 0, 2, 0),
        method = "superlevel_merge_tree",
        direction = "max"
    ))
    forest.branches <- gflow:::.basin.merge.tree.branch.table(
        forest, "max"
    )
    cross.component.ids <- vapply(
        split(
            forest.branches$basin.id,
            forest.branches$graph.component
        ),
        `[[`,
        character(1),
        1L
    )
    expect_error(
        get.basin.merge.tree.layout(
            forest,
            basin.ids = unname(cross.component.ids)
        ),
        "multiple graph components"
    )
    expect_error(
        get.basin.merge.tree.layout(
            forest,
            component = 1L,
            basin.ids = unname(cross.component.ids[[2L]])
        ),
        "do not belong"
    )
})

test_that("layout accessor rejects invalid canonical vertical values", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())

    invalid.birth <- tree
    invalid.birth$basin.table$birth.level[[1L]] <- NA_real_
    expect_error(
        get.basin.merge.tree.layout(invalid.birth),
        "branch birth, death, and persistence"
    )

    invalid.persistence <- tree
    invalid.persistence$basin.table$persistence[[1L]] <- -1
    expect_error(
        get.basin.merge.tree.layout(invalid.persistence),
        "nonnegative persistence"
    )

    invalid.event <- tree
    invalid.event$merge.table$merge.level[[1L]] <- Inf
    expect_error(
        get.basin.merge.tree.layout(invalid.event),
        "event levels and persistence"
    )
})

test_that("validated layouts reject contradictory branch and event relations", {
    for (direction in c("max", "min")) {
        tree <- get.basin.merge.tree(
            merge.tree.test.fixture(direction)
        )
        branches <- tree$basin.table[
            tree$basin.table$type == direction,
            ,
            drop = FALSE
        ]
        child.id <- branches$basin.id[
            !is.na(branches$parent.basin.id)
        ][[1L]]
        corruptions <- merge.tree.test.relational.corruptions(
            tree, direction
        )

        for (corruption.name in names(corruptions)) {
            corrupted <- corruptions[[corruption.name]]
            expect_error(
                gflow:::.validate.basin.merge.tree(corrupted),
                class = "gflow_basin_schema_error",
                info = paste(direction, corruption.name, "validator")
            )
            expect_error(
                get.basin.merge.tree.layout(
                    corrupted,
                    direction = direction
                ),
                class = "gflow_basin_schema_error",
                info = paste(direction, corruption.name, "complete")
            )
            expect_error(
                get.basin.merge.tree.layout(
                    corrupted,
                    direction = direction,
                    basin.ids = child.id,
                    close.ancestors = TRUE
                ),
                class = "gflow_basin_schema_error",
                info = paste(direction, corruption.name, "restricted")
            )
            expect_error(
                plot.basin.merge.tree(
                    corrupted,
                    direction = direction,
                    basin.ids = child.id,
                    close.ancestors = TRUE,
                    type = "tree"
                ),
                class = "gflow_basin_schema_error",
                info = paste(direction, corruption.name, "plot")
            )
        }
    }
})

test_that("relational validation preserves valid maximum and minimum layouts", {
    for (direction in c("max", "min")) {
        tree <- get.basin.merge.tree(
            merge.tree.test.fixture(direction)
        )
        complete <- get.basin.merge.tree.layout(
            tree, direction = direction
        )
        child.id <- complete$branches$basin.id[
            !is.na(complete$branches$parent.basin.id)
        ][[1L]]
        restricted <- get.basin.merge.tree.layout(
            tree,
            direction = direction,
            basin.ids = child.id,
            close.ancestors = TRUE
        )

        expect_identical(
            complete$branches,
            gflow:::.basin.merge.tree.branch.table(tree, direction)
        )
        expected.events <- tree$merge.table[
            match(
                complete$events$event.id,
                tree$merge.table$event.id
            ),
            ,
            drop = FALSE
        ]
        row.names(expected.events) <- seq_len(nrow(expected.events))
        expect_identical(complete$events, expected.events)
        expect_identical(
            restricted$branches$parent.basin.id,
            tree$basin.table$parent.basin.id[
                match(
                    restricted$branches$basin.id,
                    tree$basin.table$basin.id
                )
            ]
        )

        plot.file <- tempfile(fileext = ".pdf")
        grDevices::pdf(plot.file, width = 8, height = 6)
        plotted <- plot.basin.merge.tree(
            tree,
            direction = direction,
            basin.ids = child.id,
            close.ancestors = TRUE,
            type = "tree"
        )
        grDevices::dev.off()
        expect_identical(plotted$layout, restricted)
        expect_gt(file.info(plot.file)$size, 0)
        unlink(plot.file)
    }
})

test_that("plot consumes the public complete and restricted layout inputs", {
    tree <- get.basin.merge.tree(merge.tree.test.fixture())
    child.id <- "basin_max_v00000002"
    expect_identical(
        names(formals(plot.basin.merge.tree))[seq_len(6L)],
        c("x", "direction", "component", "type", "label", "labels")
    )
    expected <- get.basin.merge.tree.layout(
        tree,
        basin.ids = child.id,
        close.ancestors = TRUE,
        label = "extremum.vertex"
    )

    plot.file <- tempfile(fileext = ".pdf")
    grDevices::pdf(plot.file, width = 8, height = 6)
    plotted <- plot.basin.merge.tree(
        tree,
        basin.ids = child.id,
        close.ancestors = TRUE,
        type = "tree",
        label = "extremum.vertex"
    )
    grDevices::dev.off()

    expect_gt(file.info(plot.file)$size, 0)
    expect_identical(plotted$branches, expected$branches)
    expect_identical(plotted$layout$events, expected$events)
    expect_identical(
        plotted$layout$closure.added.ids,
        expected$closure.added.ids
    )
    expect_identical(plotted$layout$leaf.order, expected$leaf.order)
    expect_identical(plotted$coordinates, expected$coordinates)
    expect_identical(
        expected$closure.added.ids,
        "basin_max_v00000001"
    )
    unlink(plot.file)
})

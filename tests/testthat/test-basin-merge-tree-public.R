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

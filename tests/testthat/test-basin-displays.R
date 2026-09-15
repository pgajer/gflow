display.fixture <- function(method = "superlevel_merge_tree", direction = "max", ...) {
    a <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
    create.basin.complex(a, lapply(a, function(v) rep(1, length(v))), c(0, 3, 1, 2, 0),
                         method = method, direction = direction,
                         vertex.id = paste0("sample-", 1:5), ...)
}

test_that("summary coverage counts vertices and retains full ranking data", {
    x <- display.fixture()
    s <- summary(x, top.k.max = 1L)
    expect_equal(s$coverage$raw.vertices, 5L)
    expect_equal(s$coverage$retained.vertices, 5L)
    expect_equal(s$coverage$assigned.vertices, 5L)
    expect_equal(s$coverage$retained.basins, 2L)
    expect_equal(nrow(s$maxima), 1L)
    expect_gt(s$n.memberships, s$coverage$raw.vertices)
    text <- capture.output(print(summary(x), n = 1))
    expect_true(any(grepl("1 of 2 returned", text)))
    expect_true(any(grepl("sample-", text)))
    expect_output(print(s), "Persistence is in field units")
    expect_identical(print(s, n = 0L), s)
    expect_error(print(s, n = -1), "n")
})

test_that("raw, retained and primary coverage reflect refinement and unassigned rows", {
    x <- display.fixture(simplify.params = list(support.filter = list(enabled = TRUE, min.basin.size = 2L)))
    s <- summary(x)
    retained <- x$basin.table[x$basin.table$retained, , drop = FALSE]
    expect_equal(s$coverage$retained.vertices, length(unique(unlist(retained$retained.support.vertices))))
    expect_equal(s$coverage$assigned.vertices, sum(!is.na(x$assignment$basin.id)))
    unavailable <- display.fixture("geodesic_reachability")
    u <- summary(unavailable)
    expect_gt(u$n.assignments, 0L)
    expect_equal(u$coverage$assigned.vertices, 0L)
    expect_gt(u$coverage$raw.vertices, 0L)
})

test_that("failed objects reveal diagnostics and cannot draw analytical results", {
    x <- display.fixture()
    failed <- gflow:::.new.failed.basin.complex(
        x$method, x$direction, x$graph.input, x$field, x$parameters,
        "gflow_basin_backend_error", "Example backend unavailable")
    expect_output(print(summary(failed)), "Example backend unavailable")
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file)
    on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
    expect_warning(plot(failed, main = "Custom title"), "construction|Construction")
    expect_error(plot(failed, view = "merge_tree"), "Example backend unavailable")
    expect_error(plot(failed, view = "assignment"), "construction failed")
})

test_that("graph views preserve objects and make unavailable assignments explicit", {
    x <- display.fixture()
    coordinates <- cbind(x = 1:5, y = c(0, 1, 0, 1, 0))
    rownames(coordinates) <- x$graph.input$vertex.id
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file)
    on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
    for (view in c("field", "merge_tree", "assignment", "overlap")) {
        expect_identical(plot(x, view = view, coordinates = coordinates), x)
    }
    empty <- display.fixture(simplify.params = list(support.filter = list(enabled = TRUE, min.basin.size = 6L)))
    expect_equal(summary(empty)$coverage$retained.vertices, 0L)
    expect_identical(plot(empty, view = "overlap", coordinates = coordinates), empty)
    unavailable <- display.fixture("geodesic_reachability")
    expect_error(plot(unavailable, view = "assignment", coordinates = coordinates), "unavailable")
    expect_error(plot(x, view = "assignment"), "coordinates")
    expect_error(plot(x, view = "overlap", coordinates = coordinates[5:1, ]), "vertex IDs")
    expect_error(plot(x, view = "overlap", coordinates = coordinates, direction = "min"), "not constructed")
    both <- display.fixture(direction = "both")
    expect_error(plot(both, view = "merge_tree"), "Choose direction")
    partial <- x
    partial$status <- "partial"
    expect_warning(plot(partial), "Construction partial")
})

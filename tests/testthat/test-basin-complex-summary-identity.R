summary.identity.fixture <- function() {
    edges <- matrix(
        c(
            1, 2, 1, 1, 4, 1, 2, 3, 1, 2, 5, 1,
            3, 6, 1, 4, 5, 1, 4, 7, 1, 5, 6, 1,
            5, 8, 1, 6, 9, 1, 7, 8, 1, 8, 9, 1
        ),
        ncol = 3,
        byrow = TRUE
    )
    adjacency <- edge.lengths <- vector("list", 9)
    for (index in seq_len(nrow(edges))) {
        from <- edges[index, 1L]
        to <- edges[index, 2L]
        weight <- edges[index, 3L]
        adjacency[[from]] <- c(adjacency[[from]], to)
        adjacency[[to]] <- c(adjacency[[to]], from)
        edge.lengths[[from]] <- c(edge.lengths[[from]], weight)
        edge.lengths[[to]] <- c(edge.lengths[[to]], weight)
    }
    list(
        adj.list = adjacency,
        edge.length.list = edge.lengths,
        field = c(0, 1, 0, 1, 3, 1, 0, 1, 2),
        mass = seq_len(9),
        vertex.id = paste0("sample-", seq_len(9)),
        method.params = list(
            modulation = "CLOSEST",
            plateau.policy = "connected_exact",
            edge.length.quantile.thld = 1,
            long.edge.fallback = "allow_and_flag",
            store.trajectories = FALSE,
            symmetric.seeding = FALSE,
            tie.breaking = FALSE,
            primary.assignment.policy = "backend_primary"
        )
    )
}

summary.identity.object <- function(mass = TRUE, provenance = NULL) {
    fixture <- summary.identity.fixture()
    create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        direction = "both",
        vertex.mass = if (mass) fixture$mass else NULL,
        method.params = fixture$method.params,
        vertex.id = fixture$vertex.id,
        vertex.mass.provenance = provenance
    )
}

test_that("new constructor formals preserve complete positional compatibility", {
    expected <- c(
        "adj.list", "edge.length.list", "field", "method", "direction",
        "vertex.mass", "vertex.density", "graph.params", "method.params",
        "simplify.params", "verbose"
    )
    expect_identical(names(formals(create.basin.complex))[seq_along(expected)], expected)

    fixture <- summary.identity.fixture()
    legacy <- do.call(create.basin.complex, unname(list(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        "trajectory_flow",
        "both",
        NULL,
        NULL,
        list(edge.length.symmetry.tolerance = sqrt(.Machine$double.eps)),
        fixture$method.params,
        list(),
        FALSE
    )))
    expect_s3_class(legacy, "basin_complex")
    expect_identical(
        legacy$graph.input$vertex.id,
        as.character(seq_along(fixture$field))
    )
})

test_that("external vertex IDs are canonical and preserve internal indices", {
    object <- summary.identity.object(mass = FALSE)
    fixture <- summary.identity.fixture()
    expect_identical(
        object$assignment$vertex.id,
        fixture$vertex.id[object$assignment$vertex]
    )
    expect_identical(
        object$membership$vertex.id,
        fixture$vertex.id[object$membership$vertex]
    )
    expect_identical(
        object$basin.table$extremum.vertex.id,
        fixture$vertex.id[object$basin.table$extremum.vertex]
    )
    expect_type(object$assignment$vertex, "integer")
    expect_type(object$basin.table$raw.support.vertices, "list")

    expect_error(
        create.basin.complex(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field,
            vertex.id = factor(fixture$vertex.id)
        ),
        "must not be a factor"
    )
    for (bad in list(
        replace(fixture$vertex.id, 1L, NA_character_),
        replace(fixture$vertex.id, 1L, fixture$vertex.id[[2L]]),
        replace(fixture$vertex.id, 1L, "")
    )) {
        expect_error(
            create.basin.complex(
                fixture$adj.list,
                fixture$edge.length.list,
                fixture$field,
                vertex.id = bad
            ),
            "vertex.id"
        )
    }
    invalid <- rawToChar(as.raw(255L))
    Encoding(invalid) <- "bytes"
    expect_error(
        create.basin.complex(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field,
            vertex.id = replace(fixture$vertex.id, 1L, invalid)
        ),
        "encoding"
    )
})

test_that("mass provenance preserves verification boundaries", {
    fixture <- summary.identity.fixture()
    attestation <- list(
        claim = "occupation probability and source alignment",
        authority = "fixture manifest",
        validator = "test source contract",
        validator.version = "1",
        algorithm = "exact ordered ID and content digest comparison",
        evidence.fingerprint = "fixture-evidence-1",
        status = "validated"
    )
    provenance <- list(
        mass.kind = "occupation_probability",
        source.id = "fixture-density",
        source.fingerprint = "fixture-source-1",
        attestations = list(attestation)
    )
    object <- summary.identity.object(provenance = provenance)
    summary <- summary(object)
    stored <- summary$mass.provenance
    expect_identical(
        stored$validated.declarations$mass.kind,
        "occupation_probability"
    )
    expect_identical(
        stored$upstream.attestations[[1L]]$authority,
        "fixture manifest"
    )
    expect_true(isTRUE(stored$computed$normalized))
    expect_false("validation.status" %in% names(stored))

    provenance$declared.mass.fingerprint <- "wrong"
    expect_error(
        summary.identity.object(provenance = provenance),
        class = "gflow_basin_provenance_error"
    )

    actual <- object$provenance$mass$computed
    for (field in c(
        "declared.vertex.id.fingerprint",
        "declared.internal.graph.fingerprint"
    )) {
        mismatch <- provenance
        mismatch[[field]] <- "wrong"
        expect_error(
            summary.identity.object(provenance = mismatch),
            class = "gflow_basin_provenance_error"
        )
    }

    round.trip <- unserialize(serialize(object, NULL, version = 3L))
    expect_identical(
        round.trip$provenance$mass,
        object$provenance$mass
    )
    expect_identical(
        round.trip$graph.input$vertex.id,
        object$graph.input$vertex.id
    )
    expect_identical(
        actual$internal.graph.fingerprint,
        round.trip$provenance$mass$computed$internal.graph.fingerprint
    )
})

test_that("ranked summary resolves independently and supports zero Top-K", {
    with.mass <- summary.identity.object(mass = TRUE)
    without.mass <- summary.identity.object(mass = FALSE)
    mass.summary <- summary(with.mass, top.k.max = 0, top.k.min = 2)
    size.summary <- summary(without.mass, top.k.max = 1, top.k.min = 0)

    expect_identical(
        unname(mass.summary$rank.resolved),
        rep("primary.support.mass", 2L)
    )
    expect_identical(
        unname(size.summary$rank.resolved),
        rep("primary.support.size", 2L)
    )
    expect_equal(nrow(mass.summary$maxima), 0L)
    expect_equal(nrow(mass.summary$minima), 2L)
    expect_equal(nrow(size.summary$maxima), 1L)
    expect_equal(nrow(size.summary$minima), 0L)
    expect_true(all(
        names(mass.summary$basin.table)[
            !vapply(mass.summary$basin.table, is.list, logical(1))
        ] %in% mass.summary$column.definitions$field
    ))
})

test_that("auto distinguishes allocated and overlapping coverage mass", {
    object <- summary.identity.object(mass = TRUE)
    object$basin.table$primary.support.mass[] <- NA_real_
    object$basin.table$primary.support.size[] <- 0L
    summary <- summary(object)
    expect_identical(
        unname(summary$rank.resolved),
        rep("raw.allocated.mass", 2L)
    )

    object$provenance$allocation$raw.current <- FALSE
    stale <- summary(object)
    expect_identical(
        unname(stale$rank.resolved),
        rep("retained.support.mass", 2L)
    )
    expect_error(
        summary(object, rank.by = "raw.allocated.mass"),
        "stale_raw_allocation"
    )

    max.rows <- object$basin.table$type == "max"
    object$basin.table$retained[max.rows] <- FALSE
    empty <- summary(object)
    expect_identical(unname(empty$rank.status[["max"]]), "empty")
    expect_true(is.na(empty$rank.resolved[["max"]]))
})

test_that("mass-bearing RTCB and overlap summaries use available measures", {
    fixture <- summary.identity.fixture()
    method.params <- list(
        rtcb = list(
            edge.length.quantile.thld = 1,
            n.min = 1L,
            m.min = 0.5,
            q.min = 0.01
        ),
        overlap_cell_complex = list()
    )
    for (method in names(method.params)) {
        object <- create.basin.complex(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field,
            method = method,
            direction = "both",
            vertex.mass = fixture$mass,
            method.params = method.params[[method]],
            vertex.id = fixture$vertex.id
        )
        ranked <- summary(object)
        expect_identical(
            unname(ranked$rank.resolved),
            rep("raw.allocated.mass", 2L)
        )
        allocated.total <- vapply(
            split(object$basin.table, object$basin.table$type),
            function(rows) sum(rows$raw.allocated.mass),
            numeric(1)
        )
        expect_true(all(is.finite(allocated.total) & allocated.total > 0))
        coverage.total <- vapply(
            split(object$basin.table, object$basin.table$type),
            function(rows) sum(rows$raw.support.mass),
            numeric(1)
        )
        if (identical(method, "overlap_cell_complex")) {
            expect_equal(allocated.total, c(max = 1, min = 1))
            expect_true(any(abs(coverage.total - allocated.total) > 1e-8))
        } else {
            expect_equal(coverage.total, allocated.total)
        }
    }
})

test_that("auto rejects partial and all-zero candidates before fallback", {
    object <- summary.identity.object(mass = TRUE)
    object$basin.table$primary.support.mass[] <- 0
    max.rows <- which(object$basin.table$type == "max")
    object$basin.table$raw.allocated.mass[max.rows[[1L]]] <- NA_real_
    result <- summary(object)
    expect_identical(
        unname(result$rank.resolved[["max"]]),
        "retained.support.mass"
    )
    expect_identical(
        unname(result$rank.resolved[["min"]]),
        "raw.allocated.mass"
    )
    availability <- result$measure.availability$max
    expect_identical(
        availability$reason[
            availability$measure == "primary.support.mass"
        ],
        "all_zero"
    )
    expect_identical(
        availability$reason[
            availability$measure == "raw.allocated.mass"
        ],
        "nonfinite_or_partial"
    )
})

test_that("complete code manifests detect code changes and additions", {
    root <- tempfile("gflow-manifest-fixture-")
    dir.create(file.path(root, "R"), recursive = TRUE)
    dir.create(file.path(root, "src"), recursive = TRUE)
    writeLines("Package: gflow", file.path(root, "DESCRIPTION"))
    writeLines("exportPattern('^[^.]')", file.path(root, "NAMESPACE"))
    writeLines("f <- function() 1", file.path(root, "R", "f.R"))
    writeLines("int f() { return 1; }", file.path(root, "src", "f.cpp"))
    first <- gflow:::.gflow.compute.code.manifest(root)
    first.digest <- gflow:::.basin.hash.object(first)

    writeLines("f <- function() 2", file.path(root, "R", "f.R"))
    second <- gflow:::.gflow.compute.code.manifest(root)
    expect_false(identical(first.digest, gflow:::.basin.hash.object(second)))

    writeLines("g <- function() 3", file.path(root, "R", "g.R"))
    third <- gflow:::.gflow.compute.code.manifest(root)
    expect_true("R/g.R" %in% third$path)
    expect_false(identical(
        gflow:::.basin.hash.object(second),
        gflow:::.basin.hash.object(third)
    ))
    unlink(root, recursive = TRUE)
})

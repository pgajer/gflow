# Basin-complex identity, hashing, and provenance helpers.

.basin.hash.raw <- function(x) {
    path <- tempfile("gflow-basin-hash-")
    on.exit(unlink(path), add = TRUE)
    con <- file(path, open = "wb")
    writeBin(x, con)
    close(con)
    unname(tools::md5sum(path))
}

.basin.hash.object <- function(x) {
    .basin.hash.raw(serialize(x, NULL, version = 3L))
}

.basin.canonical.vertex.id <- function(vertex.id, n) {
    if (is.factor(vertex.id)) {
        .stop.basin.complex(
            "'vertex.id' must not be a factor; supply integer or character IDs.",
            "vertex.id"
        )
    }
    if (!(is.integer(vertex.id) || is.character(vertex.id)) ||
        is.object(vertex.id) || length(vertex.id) != n) {
        .stop.basin.complex(
            sprintf(
                "'vertex.id' must be an integer or character vector of length %d.",
                n
            ),
            "vertex.id"
        )
    }
    if (anyNA(vertex.id)) {
        .stop.basin.complex(
            "'vertex.id' must not contain missing values.",
            "vertex.id"
        )
    }
    if (is.integer(vertex.id)) {
        out <- sprintf("%d", vertex.id)
        input.type <- "integer"
    } else {
        valid <- !is.na(iconv(vertex.id, from = "", to = "UTF-8"))
        if (any(!valid)) {
            .stop.basin.complex(
                "'vertex.id' contains invalid character encoding.",
                "vertex.id"
            )
        }
        out <- enc2utf8(vertex.id)
        input.type <- "character"
    }
    if (any(!nzchar(out)) || anyDuplicated(out)) {
        .stop.basin.complex(
            "'vertex.id' must contain unique, nonempty IDs.",
            "vertex.id"
        )
    }
    list(
        values = out,
        input.type = input.type,
        encoding = "UTF-8",
        fingerprint = .basin.hash.object(out)
    )
}

.basin.internal.graph.fingerprint <- function(graph, vertex.id) {
    .basin.hash.object(list(
        vertex.id = vertex.id,
        adj.list = lapply(graph$adj.list, as.integer),
        edge.length.list = lapply(graph$edge.length.list, as.numeric)
    ))
}

.basin.external.id.values <- function(index, vertex.id) {
    index <- suppressWarnings(as.integer(index))
    out <- rep.int(NA_character_, length(index))
    valid <- !is.na(index) & index >= 1L & index <= length(vertex.id)
    out[valid] <- vertex.id[index[valid]]
    out
}

.basin.external.id.lists <- function(indices, vertex.id) {
    lapply(indices, function(index) {
        .basin.external.id.values(index, vertex.id)
    })
}

.basin.attach.external.vertex.ids <- function(extrema,
                                              basin.table,
                                              membership,
                                              assignment,
                                              merge.table,
                                              vertex.id) {
    extrema$representative.vertex.id <- .basin.external.id.values(
        extrema$representative.vertex,
        vertex.id
    )
    extrema$plateau.vertex.ids <- I(.basin.external.id.lists(
        extrema$plateau.vertices,
        vertex.id
    ))

    basin.table$extremum.vertex.id <- .basin.external.id.values(
        basin.table$extremum.vertex,
        vertex.id
    )
    for (prefix in c("raw.support", "retained.support", "primary.support")) {
        internal <- paste0(prefix, ".vertices")
        external <- paste0(prefix, ".vertex.ids")
        basin.table[[external]] <- I(.basin.external.id.lists(
            basin.table[[internal]],
            vertex.id
        ))
    }

    membership$vertex.id <- .basin.external.id.values(
        membership$vertex,
        vertex.id
    )
    assignment$vertex.id <- .basin.external.id.values(
        assignment$vertex,
        vertex.id
    )
    assignment$root.vertex.id <- .basin.external.id.values(
        assignment$root.vertex,
        vertex.id
    )
    assignment$next.vertex.id <- .basin.external.id.values(
        assignment$next.vertex,
        vertex.id
    )

    if (!is.null(merge.table)) {
        merge.table$merge.vertex.ids <- I(.basin.external.id.lists(
            merge.table$merge.vertices,
            vertex.id
        ))
    }
    list(
        extrema = extrema,
        basin.table = basin.table,
        membership = membership,
        assignment = assignment,
        merge.table = merge.table
    )
}

.basin.mass.kinds <- c(
    "occupation_probability",
    "empirical_weight",
    "quadrature_weight",
    "population_weight",
    "unspecified_explicit"
)

.basin.attestation.statuses <- c("validated", "not_checked")

.basin.validate.attestation.evidence <- function(x, argument) {
    if (is.null(x)) {
        return(NULL)
    }
    x <- .basin.assert.named.list(x, argument)
    if (anyDuplicated(names(x)) || any(!nzchar(names(x)))) {
        .stop.basin.complex(
            sprintf("'%s' must have unique, nonempty field names.", argument),
            argument
        )
    }
    supported <- vapply(x, function(value) {
        (is.character(value) || is.integer(value) ||
            is.numeric(value) || is.logical(value)) &&
            length(value) == 1L && !is.na(value) &&
            if (is.character(value)) nzchar(trimws(value)) else is.finite(value)
    }, logical(1))
    if (any(!supported)) {
        .stop.basin.complex(
            sprintf(
                "'%s' values must be finite, nonmissing atomic scalars.",
                argument
            ),
            argument
        )
    }
    x[] <- lapply(x, function(value) {
        if (is.character(value)) enc2utf8(value) else value
    })
    x
}

.basin.validate.attestation <- function(x, index) {
    argument <- sprintf("vertex.mass.provenance$attestations[[%d]]", index)
    x <- .basin.assert.named.list(x, argument)
    required <- c(
        "claim", "authority", "validator", "validator.version",
        "algorithm", "evidence.fingerprint", "status"
    )
    missing <- setdiff(required, names(x))
    if (length(missing) > 0L) {
        .stop.basin.complex(
            sprintf(
                "'%s' is missing: %s.",
                argument,
                paste(missing, collapse = ", ")
            ),
            argument
        )
    }
    out <- lapply(required, function(name) {
        value <- x[[name]]
        if (!is.character(value) || length(value) != 1L || is.na(value) ||
            !nzchar(trimws(value))) {
            .stop.basin.complex(
                sprintf("'%s$%s' must be one nonempty string.", argument, name),
                paste0(argument, "$", name)
            )
        }
        enc2utf8(value)
    })
    names(out) <- required
    out$status <- .basin.assert.choice(
        out$status,
        .basin.attestation.statuses,
        paste0(argument, "$status")
    )
    optional.strings <- c(
        "contract.version", "source.id", "source.fingerprint"
    )
    for (name in optional.strings) {
        value <- x[[name]]
        if (is.null(value)) {
            next
        }
        if (!is.character(value) || length(value) != 1L || is.na(value) ||
            !nzchar(trimws(value))) {
            .stop.basin.complex(
                sprintf("'%s$%s' must be one nonempty string.", argument, name),
                paste0(argument, "$", name)
            )
        }
        out[[name]] <- enc2utf8(value)
    }
    out$evidence <- .basin.validate.attestation.evidence(
        x$evidence,
        paste0(argument, "$evidence")
    )
    out
}

.basin.validate.mass.provenance <- function(provenance,
                                            vertex.mass,
                                            mass.total,
                                            normalized.mass,
                                            vertex.identity,
                                            internal.graph.fingerprint) {
    if (is.null(vertex.mass)) {
        if (!is.null(provenance)) {
            .stop.basin.complex(
                "'vertex.mass.provenance' must be NULL when 'vertex.mass' is NULL.",
                "vertex.mass.provenance"
            )
        }
        return(NULL)
    }

    supplied <- !is.null(provenance)
    if (!supplied) {
        provenance <- list(
            mass.kind = "unspecified_explicit",
            attestations = list()
        )
    }
    provenance <- .basin.assert.named.list(
        provenance,
        "vertex.mass.provenance"
    )
    allowed <- c(
        "mass.kind",
        "source.id",
        "source.fingerprint",
        "declared.mass.fingerprint",
        "declared.vertex.id.fingerprint",
        "declared.internal.graph.fingerprint",
        "attestations"
    )
    unknown <- setdiff(names(provenance), allowed)
    if (length(unknown) > 0L) {
        .stop.basin.complex(
            sprintf(
                "Unknown vertex.mass.provenance field%s: %s.",
                if (length(unknown) == 1L) "" else "s",
                paste(unknown, collapse = ", ")
            ),
            "vertex.mass.provenance"
        )
    }

    mass.kind <- .basin.assert.choice(
        provenance$mass.kind %||% "unspecified_explicit",
        .basin.mass.kinds,
        "vertex.mass.provenance$mass.kind"
    )
    scalar.string <- function(name, required = FALSE) {
        value <- provenance[[name]]
        if (is.null(value) && !required) {
            return(NULL)
        }
        if (!is.character(value) || length(value) != 1L || is.na(value) ||
            !nzchar(trimws(value))) {
            .stop.basin.complex(
                sprintf(
                    "'vertex.mass.provenance$%s' must be one nonempty string.",
                    name
                ),
                paste0("vertex.mass.provenance$", name)
            )
        }
        enc2utf8(value)
    }

    actual.mass.fingerprint <- .basin.hash.object(as.numeric(vertex.mass))
    declared.mass <- scalar.string("declared.mass.fingerprint")
    declared.vertex <- scalar.string("declared.vertex.id.fingerprint")
    declared.graph <- scalar.string("declared.internal.graph.fingerprint")
    mismatches <- c(
        if (!is.null(declared.mass) &&
            !identical(declared.mass, actual.mass.fingerprint)) {
            "declared.mass.fingerprint"
        },
        if (!is.null(declared.vertex) &&
            !identical(declared.vertex, vertex.identity$fingerprint)) {
            "declared.vertex.id.fingerprint"
        },
        if (!is.null(declared.graph) &&
            !identical(declared.graph, internal.graph.fingerprint)) {
            "declared.internal.graph.fingerprint"
        }
    )
    if (length(mismatches) > 0L) {
        .stop.basin.complex(
            sprintf(
                "Recomputable provenance mismatch: %s.",
                paste(mismatches, collapse = ", ")
            ),
            "vertex.mass.provenance",
            class = "gflow_basin_provenance_error",
            details = list(mismatches = mismatches)
        )
    }

    attestations <- provenance$attestations %||% list()
    if (!is.list(attestations) || is.object(attestations)) {
        .stop.basin.complex(
            "'vertex.mass.provenance$attestations' must be an ordinary list.",
            "vertex.mass.provenance$attestations"
        )
    }
    attestations <- lapply(
        seq_along(attestations),
        function(index) .basin.validate.attestation(attestations[[index]], index)
    )

    list(
        computed = list(
            mass.fingerprint = actual.mass.fingerprint,
            input.total = as.numeric(mass.total),
            normalized.total = as.numeric(sum(normalized.mass)),
            normalized = isTRUE(all.equal(
                sum(normalized.mass),
                1,
                tolerance = 1e-12
            )),
            vertex.id.fingerprint = vertex.identity$fingerprint,
            internal.graph.fingerprint = internal.graph.fingerprint
        ),
        validated.declarations = list(
            mass.kind = mass.kind,
            declared.mass.fingerprint = declared.mass,
            declared.vertex.id.fingerprint = declared.vertex,
            declared.internal.graph.fingerprint = declared.graph,
            schema = "gflow_basin_mass_provenance/1"
        ),
        upstream.declarations = list(
            source.id = scalar.string("source.id"),
            source.fingerprint = scalar.string("source.fingerprint"),
            status = "not_checked"
        ),
        upstream.attestations = attestations,
        supplied = supplied
    )
}

.gflow.code.manifest.version <- "1"
.gflow.build.identity.cache <- new.env(parent = emptyenv())

.gflow.package.root <- function() {
    path <- tryCatch(
        getNamespaceInfo(asNamespace("gflow"), "path"),
        error = function(e) ""
    )
    path <- as.character(path %||% "")
    if (length(path) == 1L && nzchar(path) &&
        file.exists(file.path(path, "DESCRIPTION"))) {
        return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
    ""
}

.gflow.code.input.files <- function(root) {
    roots <- c("R", "src", file.path("inst", "include"))
    files <- unlist(lapply(roots, function(relative) {
        path <- file.path(root, relative)
        if (!dir.exists(path)) {
            return(character())
        }
        list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    }), use.names = FALSE)
    top <- c(
        "DESCRIPTION", "NAMESPACE", "configure", "configure.ac",
        "cleanup", "cleanup.win"
    )
    files <- c(files, file.path(root, top[file.exists(file.path(root, top))]))
    files <- files[file.info(files)$isdir %in% FALSE]
    files <- files[!grepl(
        "(^|/)(gflow-code-manifest\\.tsv|\\.DS_Store)$|\\.(o|so|dll|dylib)$",
        files
    )]
    sort(unique(normalizePath(files, winslash = "/", mustWork = TRUE)))
}

.gflow.compute.code.manifest <- function(root) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    files <- .gflow.code.input.files(root)
    relative <- substring(files, nchar(root) + 2L)
    data.frame(
        path = relative,
        md5 = unname(tools::md5sum(files)),
        stringsAsFactors = FALSE
    )
}

.gflow.embedded.manifest.path <- function(root) {
    file.path(root, "extdata", "gflow-code-manifest.tsv")
}

.gflow.write.code.manifest <- function(root = ".") {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    manifest <- .gflow.compute.code.manifest(root)
    path <- file.path(root, "inst", "extdata", "gflow-code-manifest.tsv")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
        manifest,
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        fileEncoding = "UTF-8"
    )
    invisible(path)
}

.gflow.read.code.manifest <- function(root) {
    if (dir.exists(file.path(root, "R"))) {
        return(.gflow.compute.code.manifest(root))
    }
    path <- .gflow.embedded.manifest.path(root)
    if (!file.exists(path)) {
        return(data.frame(path = character(), md5 = character()))
    }
    utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        encoding = "UTF-8"
    )
}

.gflow.native.artifact.digest <- function(root) {
    files <- list.files(
        root,
        pattern = "\\.(so|dll|dylib)$",
        recursive = TRUE,
        full.names = TRUE
    )
    files <- sort(files[file.info(files)$isdir %in% FALSE])
    if (length(files) == 0L) {
        return(NA_character_)
    }
    .basin.hash.object(stats::setNames(
        unname(tools::md5sum(files)),
        substring(files, nchar(root) + 2L)
    ))
}

.gflow.git.value <- function(root, args) {
    git.marker <- file.path(root, ".git")
    if (!file.exists(git.marker) && !dir.exists(git.marker)) {
        return(character())
    }
    tryCatch(
        trimws(suppressWarnings(system2(
            "git",
            c("-C", shQuote(root), args),
            stdout = TRUE,
            stderr = FALSE
        ))),
        error = function(e) character()
    )
}

.gflow.runtime.provenance <- function(root) {
    description <- tryCatch(
        if (file.exists(file.path(root, "DESCRIPTION"))) {
            read.dcf(file.path(root, "DESCRIPTION"))[1L, , drop = TRUE]
        } else {
            unlist(utils::packageDescription("gflow"), use.names = TRUE)
        },
        error = function(e) character()
    )
    dependency.fields <- intersect(c("Imports", "LinkingTo"), names(description))
    dependency.text <- paste(description[dependency.fields], collapse = ",")
    dependencies <- trimws(unlist(strsplit(dependency.text, ",", fixed = TRUE)))
    dependencies <- sub("\\s*\\(.*$", "", dependencies)
    dependencies <- sort(unique(dependencies[nzchar(dependencies)]))
    dependency.versions <- vapply(dependencies, function(package) {
        tryCatch(
            as.character(utils::packageVersion(package)),
            error = function(e) NA_character_
        )
    }, character(1))
    list(
        r.version = paste(R.version$major, R.version$minor, sep = "."),
        platform = R.version$platform,
        architecture = R.version$arch,
        os = R.version$os,
        dependency.versions = dependency.versions
    )
}

#' Report the Loaded gflow Build Identity
#'
#' Returns immutable package-code and native-artifact identity plus the
#' conservative runtime compatibility identity used by basin-complex caches.
#'
#' @param refresh Recompute the identity instead of returning the identity
#'   cached when the currently loaded namespace first constructed a basin.
#' @return A named list containing build and runtime identity records.
#' @export
get.gflow.build.identity <- function(refresh = FALSE) {
    refresh <- .basin.assert.logical(refresh, "refresh")
    if (!refresh && exists(
        "value",
        envir = .gflow.build.identity.cache,
        inherits = FALSE
    )) {
        return(get(
            "value",
            envir = .gflow.build.identity.cache,
            inherits = FALSE
        ))
    }
    root <- .gflow.package.root()
    if (!nzchar(root)) {
        root <- normalizePath(".", winslash = "/", mustWork = TRUE)
    }
    manifest <- .gflow.read.code.manifest(root)
    manifest.digest <- .basin.hash.object(list(
        version = .gflow.code.manifest.version,
        path = as.character(manifest$path),
        md5 = as.character(manifest$md5)
    ))
    revision <- .gflow.git.value(root, c("rev-parse", "HEAD"))
    revision <- if (length(revision) == 1L && nzchar(revision)) {
        revision
    } else {
        NA_character_
    }
    code.paths <- as.character(manifest$path)
    dirty <- if (!is.na(revision) && length(code.paths) > 0L) {
        length(.gflow.git.value(
            root,
            c("status", "--porcelain", "--", shQuote(code.paths))
        )) > 0L
    } else {
        NA
    }
    native.digest <- .gflow.native.artifact.digest(root)
    package.version <- .basin.package.version()
    build.id <- .basin.hash.object(list(
        manifest.version = .gflow.code.manifest.version,
        manifest.digest = manifest.digest,
        native.artifact.digest = native.digest,
        source.revision = revision,
        source.code.dirty = dirty,
        package.version = package.version
    ))
    runtime <- .gflow.runtime.provenance(root)
    runtime$id <- .basin.hash.object(runtime)
    result <- list(
        package.version = package.version,
        manifest.version = .gflow.code.manifest.version,
        manifest.digest = manifest.digest,
        native.artifact.digest = native.digest,
        source.revision = revision,
        source.code.dirty = dirty,
        build.id = build.id,
        runtime = runtime
    )
    assign("value", result, envir = .gflow.build.identity.cache)
    result
}

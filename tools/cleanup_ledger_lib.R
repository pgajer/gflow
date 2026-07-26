cleanup.repo.root <- function(start = getwd()) {
    path <- normalizePath(start, mustWork = TRUE)
    repeat {
        if (file.exists(file.path(path, "DESCRIPTION")) &&
            file.exists(file.path(path, "NAMESPACE"))) {
            return(path)
        }
        parent <- dirname(path)
        if (identical(parent, path)) {
            stop("Cannot locate the gflow package root.", call. = FALSE)
        }
        path <- parent
    }
}

cleanup.namespace.inventory <- function(root) {
    lines <- trimws(readLines(file.path(root, "NAMESPACE"), warn = FALSE))
    export.lines <- grep("^export\\(", lines, value = TRUE)
    s3.lines <- grep("^S3method\\(", lines, value = TRUE)

    exports <- substring(export.lines, 8L, nchar(export.lines) - 1L)
    s3.payload <- substring(s3.lines, 10L, nchar(s3.lines) - 1L)
    s3.parts <- strsplit(s3.payload, ",", fixed = TRUE)
    s3.generic <- vapply(s3.parts, function(x) trimws(x[[1L]]), character(1))
    s3.class <- vapply(s3.parts, function(x) trimws(x[[2L]]), character(1))

    rbind(
        data.frame(
            item_type = "export",
            item_name = exports,
            generic = "",
            class = "",
            key = paste0("export:", exports),
            stringsAsFactors = FALSE
        ),
        data.frame(
            item_type = "s3",
            item_name = paste0(s3.generic, ".", s3.class),
            generic = s3.generic,
            class = s3.class,
            key = paste0("s3:", s3.generic, ":", s3.class),
            stringsAsFactors = FALSE
        )
    )
}

cleanup.function.definitions <- function(root) {
    files <- sort(list.files(file.path(root, "R"), "\\.[Rr]$",
                             full.names = TRUE))
    definitions <- list()

    for (file in files) {
        expressions <- parse(file, keep.source = TRUE)
        for (expression in expressions) {
            is.assignment <- is.call(expression) &&
                (identical(expression[[1L]], as.name("<-")) ||
                 identical(expression[[1L]], as.name("=")))
            if (!is.assignment || !is.symbol(expression[[2L]]) ||
                !is.call(expression[[3L]]) ||
                !identical(expression[[3L]][[1L]], as.name("function"))) {
                next
            }
            name <- as.character(expression[[2L]])
            definitions[[name]] <- list(
                name = name,
                file = sub(paste0("^", root, "/"), "", file),
                expression = expression[[3L]]
            )
        }
    }
    definitions
}

cleanup.call.symbols <- function(expression) {
    calls <- character()
    walk <- function(node) {
        if (!is.call(node)) {
            return(invisible(NULL))
        }
        head <- node[[1L]]
        if (is.symbol(head)) {
            calls <<- c(calls, as.character(head))
        }
        for (i in seq_along(node)[-1L]) {
            walk(node[[i]])
        }
        invisible(NULL)
    }
    walk(expression)
    ignored <- c(
        "{", "(", "[", "[[", "$", "@", "<-", "=", "<<-", "~",
        "if", "for", "while", "repeat", "break", "next", "return",
        "function", "::", ":::", "+", "-", "*", "/", "^", "%%", "%/%",
        "%*%", "%in%", ":", "==", "!=", "<", ">", "<=", ">=", "&", "&&",
        "|", "||", "!", "?"
    )
    sort(unique(setdiff(calls, ignored)))
}

cleanup.call.graph <- function(definitions) {
    names.all <- names(definitions)
    forward <- setNames(vector("list", length(names.all)), names.all)
    for (name in names.all) {
        calls <- cleanup.call.symbols(definitions[[name]]$expression)
        forward[[name]] <- sort(intersect(calls, names.all))
    }
    reverse <- setNames(vector("list", length(names.all)), names.all)
    for (caller in names.all) {
        for (callee in forward[[caller]]) {
            reverse[[callee]] <- c(reverse[[callee]], caller)
        }
    }
    reverse <- lapply(reverse, function(x) sort(unique(x)))
    list(forward = forward, reverse = reverse)
}

cleanup.transitive.callers <- function(name, reverse) {
    pending <- reverse[[name]]
    found <- character()
    while (length(pending)) {
        current <- pending[[1L]]
        pending <- pending[-1L]
        if (current %in% found) {
            next
        }
        found <- c(found, current)
        pending <- c(pending, reverse[[current]])
    }
    sort(unique(found))
}

cleanup.transitive.dependencies <- function(roots, forward) {
    pending <- intersect(roots, names(forward))
    found <- character()
    while (length(pending)) {
        current <- pending[[1L]]
        pending <- pending[-1L]
        if (current %in% found) {
            next
        }
        found <- c(found, current)
        pending <- c(pending, forward[[current]])
    }
    sort(unique(found))
}

cleanup.primary.protected.files <- function() {
    c(
        "R/2d_smooth_morse_smale.R",
        "R/basin_complex.R",
        "R/basin_complex_converters.R",
        "R/basin_cx.R",
        "R/cluster_local_extrema.R",
        "R/compute_basins_of_attraction.R",
        "R/compute_gfc.R",
        "R/compute_gfc_basins.R",
        "R/compute_partition_cell_stats.R",
        "R/compute_pextrema_nbhds.R",
        "R/extrema_utils.R",
        "R/filter_basins.R",
        "R/filter_basins_by_relvalue.R",
        "R/gfc_flow.R",
        "R/gflow_basins.R",
        "R/gflow_cx.R",
        "R/gradient_basin_utilities.R",
        "R/graph_gradient_flow.R",
        "R/graph_gradient_flow_basins.R",
        "R/graph_MS_cx.R",
        "R/local_extrema.R",
        "R/merge_clustered_basins.R",
        "R/pextrema.R"
    )
}

cleanup.protected.roots <- function(definitions) {
    files <- cleanup.primary.protected.files()
    names(Filter(
        function(definition) definition$file %in% files,
        definitions
    ))
}

cleanup.protected.symbols <- function(definitions, graph) {
    cleanup.transitive.dependencies(
        cleanup.protected.roots(definitions),
        graph$forward
    )
}

cleanup.classification <- function(name, source.file, protected) {
    if (name %in% protected) {
        return(c(
            ownership = "PROTECTED",
            intended_package = "gflow",
            lifecycle_action = "PROTECT",
            target_release = "outside-cleanup-scope"
        ))
    }

    association.files <- c(
        "R/gfassoc_utils.R", "R/gfcor.R", "R/lcor.R",
        "R/lcor_with_posterior.R", "R/lslope.R",
        "R/permutation_test_lcor.R", "R/pvalue_utils.R"
    )
    if (source.file %in% association.files) {
        return(c(
            ownership = "CORE-ASSOCIATION",
            intended_package = "gflow",
            lifecycle_action = "KEEP-AND-CONSOLIDATE",
            target_release = "cleanup-release"
        ))
    }

    if (source.file == "R/compositional_data_utils.R") {
        return(c(
            ownership = "DOMAIN",
            intended_package = "microbiome.utils",
            lifecycle_action = "DEEXPORT-THEN-RELOCATE",
            target_release = "phase-5"
        ))
    }
    if (source.file == "R/two_factor_analysis.R") {
        owner <- if (grepl("cst", name, fixed = TRUE)) "gcstflow" else
            "external-analysis-repository"
        return(c(
            ownership = "DOMAIN",
            intended_package = owner,
            lifecycle_action = "DEEXPORT-THEN-RELOCATE",
            target_release = "phase-5"
        ))
    }
    if (source.file %in% c(
        "R/preprocess_metabolon_twoplane.R",
        "R/visualize_mgcp_cluster_single.R",
        "R/visualize_mgcp_clusters.R",
        "R/plot_feature_distance_diagnostics.R",
        "R/partition_heatmap.R"
    )) {
        owner <- if (grepl("visualize_mgcp|partition_heatmap", source.file))
            "gflowui" else "external-analysis-repository"
        return(c(
            ownership = "DOMAIN",
            intended_package = owner,
            lifecycle_action = "DEEXPORT-THEN-RELOCATE",
            target_release = "phase-5"
        ))
    }

    if (source.file %in% c(
        "R/diffusion_pseudotime_sparse.R", "R/phate_core.R",
        "R/potential_metric_graph.R", "R/trajectory_distances.R",
        "R/find_shortest_paths_within_radius.R", "R/graph_utils.R",
        "R/graph_endpoint_geometry.R", "R/graphics.R",
        "R/partition_graph.R"
    )) {
        temporary.public <- name %in% c(
            "phate", "phate.core", "phate.embed",
            "compute.diffusion.pseudotime.sparse",
            "compute.potential.pseudotime.sparse"
        )
        return(c(
            ownership = "DGRAPHS",
            intended_package = "dgraphs",
            lifecycle_action = if (temporary.public)
                "TEMPORARY-PUBLIC-MIGRATE" else
                "REVIEW-DEEXPORT-OR-ADAPTER",
            target_release = if (temporary.public)
                "after-dgraphs-parity" else
                "phase-5-or-6"
        ))
    }

    if (source.file == "R/clustering.R" && name %in% c(
        "cluster.graph.louvain", "congruence.with.labels"
    )) {
        return(c(
            ownership = "DGRAPHS",
            intended_package = "dgraphs",
            lifecycle_action = "TEMPORARY-PUBLIC-MIGRATE",
            target_release = "after-dgraphs-parity"
        ))
    }

    if (source.file %in% c("R/graph_generators.R", "R/hHN_graphs.R")) {
        return(c(
            ownership = "EXAMPLE",
            intended_package = "tests-or-external-examples",
            lifecycle_action = "DEEXPORT-THEN-MOVE",
            target_release = "phase-5"
        ))
    }

    if (source.file %in% c(
        "R/annotation.R", "R/circular_data_utils.R", "R/clustering.R",
        "R/density.R", "R/select_pts.R", "R/utils.R"
    )) {
        return(c(
            ownership = "CORE-PRIVATE",
            intended_package = "gflow-internal-or-specialist-package",
            lifecycle_action = "REVIEW-AND-DEEXPORT",
            target_release = "phase-5-or-6"
        ))
    }

    c(
        ownership = "CORE-ANALYSIS",
        intended_package = "gflow",
        lifecycle_action = "KEEP-PENDING-S3-REVIEW",
        target_release = "phase-6"
    )
}

cleanup.find.files <- function(root, directory, pattern) {
    path <- file.path(root, directory)
    if (!dir.exists(path)) {
        return(character())
    }
    files <- list.files(path, pattern, full.names = TRUE, recursive = TRUE)
    sub(paste0("^", root, "/"), "", files)
}

cleanup.files.containing <- function(root, files, token) {
    hits <- character()
    for (file in files) {
        lines <- readLines(file.path(root, file), warn = FALSE)
        if (any(grepl(token, lines, fixed = TRUE))) {
            hits <- c(hits, file)
        }
    }
    sort(unique(hits))
}

cleanup.documentation <- function(root, name) {
    files <- cleanup.find.files(root, "man", "\\.Rd$")
    alias <- paste0("\\alias{", name, "}")
    hits <- cleanup.files.containing(root, files, alias)
    if (length(hits)) paste(hits, collapse = ";") else "none"
}

cleanup.tests <- function(root, name) {
    files <- cleanup.find.files(root, "tests/testthat", "\\.[Rr]$")
    hits <- cleanup.files.containing(root, files, name)
    if (length(hits)) paste(hits, collapse = ";") else "none"
}

cleanup.read.downstream.registry <- function(root) {
    path <- file.path(
        root, "split_audit/cleanup/downstream-repositories.txt"
    )
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(trimws(lines)) & !grepl("^#", trimws(lines))]
    parts <- strsplit(lines, "\t", fixed = TRUE)
    data.frame(
        path = vapply(parts, `[[`, character(1), 1L),
        scope = vapply(parts, `[[`, character(1), 2L),
        reason = vapply(parts, `[[`, character(1), 3L),
        stringsAsFactors = FALSE
    )
}

cleanup.downstream.calls <- function(root, registry) {
    result <- list()
    repositories <- registry$path[
        dir.exists(registry$path) &
            normalizePath(registry$path) != normalizePath(root)
    ]
    if (!length(repositories)) {
        return(result)
    }
    if (!nzchar(Sys.which("rg"))) {
        stop("The downstream snapshot requires ripgrep (`rg`).",
             call. = FALSE)
    }

    expression <- "gflow:::{0,1}[A-Za-z][A-Za-z0-9._]*"
    args <- c(
        "-n", "-o", "--no-heading", "--hidden",
        "-g", shQuote("*.R"),
        "-g", shQuote("*.Rmd"),
        "-g", shQuote("*.qmd"),
        "-g", shQuote("*.Rnw"),
        "-g", shQuote("*.md"),
        "-g", shQuote("!**/.git/**"),
        "-g", shQuote("!**/node_modules/**"),
        "-g", shQuote("!**/gflow.Rcheck/**"),
        "-g", shQuote("!**/site_libs/**"),
        shQuote(expression),
        shQuote(repositories)
    )
    output <- suppressWarnings(system2(
        "rg", args, stdout = TRUE, stderr = FALSE
    ))
    if (!length(output)) {
        return(result)
    }

    token <- sub("^.*:[0-9]+:", "", output)
    file.line <- sub(":(gflow:::{0,1}[A-Za-z][A-Za-z0-9._]*)$", "",
                     output)
    line.number <- sub("^.*:([0-9]+)$", "\\1", file.line)
    file <- sub(":[0-9]+$", "", file.line)
    functions <- sub("^gflow:::{0,1}", "", token)

    for (i in seq_along(output)) {
        repository <- repositories[
            startsWith(file[[i]], paste0(repositories, "/"))
        ][[1L]]
        relative <- substring(file[[i]], nchar(repository) + 2L)
        reference <- paste0(
            basename(repository), ":", relative, ":", line.number[[i]]
        )
        name <- functions[[i]]
        result[[name]] <- c(result[[name]], reference)
    }
    lapply(result, function(x) sort(unique(x)))
}

cleanup.collapse <- function(x) {
    if (!length(x)) "none" else paste(sort(unique(x)), collapse = ";")
}

cleanup.build.ledger <- function(root, scan.downstream = TRUE) {
    inventory <- cleanup.namespace.inventory(root)
    definitions <- cleanup.function.definitions(root)
    graph <- cleanup.call.graph(definitions)
    protected <- cleanup.protected.symbols(definitions, graph)
    downstream <- if (scan.downstream) {
        cleanup.downstream.calls(
            root, cleanup.read.downstream.registry(root)
        )
    } else {
        list()
    }

    rows <- lapply(seq_len(nrow(inventory)), function(i) {
        item <- inventory[i, ]
        method.name <- item$item_name
        definition <- definitions[[method.name]]
        source.file <- if (is.null(definition)) "unresolved" else
            definition$file
        classification <- cleanup.classification(
            method.name, source.file, protected
        )
        direct <- graph$reverse[[method.name]]
        transitive <- cleanup.transitive.callers(method.name, graph$reverse)
        known.downstream <- downstream[[method.name]]

        data.frame(
            key = item$key,
            item_type = item$item_type,
            item_name = item$item_name,
            generic = item$generic,
            class = item$class,
            source_file = source.file,
            ownership = unname(classification[["ownership"]]),
            intended_package =
                unname(classification[["intended_package"]]),
            current_status = if (item$item_type == "export")
                "public-export" else "registered-s3",
            direct_callers = cleanup.collapse(direct),
            transitive_callers = cleanup.collapse(transitive),
            tests = cleanup.tests(root, method.name),
            documentation = cleanup.documentation(root, method.name),
            known_downstream_callers =
                cleanup.collapse(known.downstream),
            lifecycle_action =
                unname(classification[["lifecycle_action"]]),
            target_release =
                unname(classification[["target_release"]]),
            protected_call_graph =
                if (method.name %in% protected) "yes" else "no",
            stringsAsFactors = FALSE
        )
    })
    ledger <- do.call(rbind, rows)
    ledger[order(ledger$item_type, ledger$item_name), ]
}

cleanup.hash.expression <- function(expression) {
    path <- tempfile("gflow-protected-expression-", fileext = ".R")
    on.exit(unlink(path), add = TRUE)
    writeLines(
        deparse(expression, control = c("keepInteger", "keepNA")),
        path,
        useBytes = TRUE
    )
    unname(tools::md5sum(path))
}

cleanup.native.symbols <- function(protected, definitions) {
    symbols <- character()
    for (name in intersect(protected, names(definitions))) {
        text <- paste(deparse(definitions[[name]]$expression), collapse = "\n")
        hits <- regmatches(
            text,
            gregexpr("_gflow_[A-Za-z0-9_]+", text, perl = TRUE)
        )[[1L]]
        if (length(hits) && !identical(hits, character(0))) {
            symbols <- c(symbols, hits)
        }
    }
    sort(unique(symbols))
}

cleanup.build.protected.surface <- function(root) {
    definitions <- cleanup.function.definitions(root)
    graph <- cleanup.call.graph(definitions)
    protected <- cleanup.protected.symbols(definitions, graph)
    primary.files <- cleanup.primary.protected.files()
    inventory <- cleanup.namespace.inventory(root)

    lines <- c(
        "# gflow protected basin/extrema construction surface",
        "# Format: TYPE<TAB>identifier<TAB>source-or-class<TAB>md5-when-applicable",
        "# Generated by tools/build_cleanup_ledger.R; update only during explicit basin work."
    )
    for (file in primary.files) {
        hash <- unname(tools::md5sum(file.path(root, file)))
        lines <- c(lines, paste("FILE", file, "-", hash, sep = "\t"))
    }

    mixed.symbols <- protected[vapply(
        protected,
        function(name) {
            definition <- definitions[[name]]
            !is.null(definition) && !(definition$file %in% primary.files)
        },
        logical(1)
    )]
    for (name in mixed.symbols) {
        definition <- definitions[[name]]
        lines <- c(lines, paste(
            "SYMBOL", name, definition$file,
            cleanup.hash.expression(definition$expression),
            sep = "\t"
        ))
    }

    protected.exports <- inventory$item_name[
        inventory$item_type == "export" &
            inventory$item_name %in% protected
    ]
    for (name in sort(protected.exports)) {
        lines <- c(lines, paste("EXPORT", name, sep = "\t"))
    }

    protected.s3 <- inventory[
        inventory$item_type == "s3" &
            inventory$item_name %in% protected,
    ]
    if (nrow(protected.s3)) {
        for (i in seq_len(nrow(protected.s3))) {
            lines <- c(lines, paste(
                "S3", protected.s3$generic[[i]],
                protected.s3$class[[i]], sep = "\t"
            ))
        }
    }

    for (symbol in cleanup.native.symbols(protected, definitions)) {
        lines <- c(lines, paste(
            "NATIVE", symbol, "R/RcppExports.R", sep = "\t"
        ))
    }
    lines
}

cleanup.read.protected.surface <- function(root) {
    path <- file.path(
        root, "split_audit/cleanup/protected-basin-surface.txt"
    )
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(trimws(lines)) & !grepl("^#", trimws(lines))]
    parts <- strsplit(lines, "\t", fixed = TRUE)
    width <- max(lengths(parts))
    parts <- lapply(parts, function(x) c(x, rep("", width - length(x))))
    matrix <- do.call(rbind, parts)
    data.frame(
        type = matrix[, 1L],
        identifier = matrix[, 2L],
        source_or_class = matrix[, 3L],
        md5 = matrix[, 4L],
        stringsAsFactors = FALSE
    )
}

cleanup.check.protected.surface <- function(root) {
    surface <- cleanup.read.protected.surface(root)
    definitions <- cleanup.function.definitions(root)
    errors <- character()

    for (i in seq_len(nrow(surface))) {
        row <- surface[i, ]
        if (row$type == "FILE") {
            path <- file.path(root, row$identifier)
            actual <- if (file.exists(path))
                unname(tools::md5sum(path)) else "<missing>"
            if (!identical(actual, row$md5)) {
                errors <- c(errors, paste("Protected file changed:",
                                          row$identifier))
            }
        } else if (row$type == "SYMBOL") {
            definition <- definitions[[row$identifier]]
            actual <- if (is.null(definition)) "<missing>" else
                cleanup.hash.expression(definition$expression)
            if (!identical(actual, row$md5)) {
                errors <- c(errors, paste("Protected symbol changed:",
                                          row$identifier))
            }
        }
    }
    errors
}

cleanup.validate.ledger <- function(root, ledger = NULL) {
    if (is.null(ledger)) {
        ledger <- read.csv(
            file.path(root, "split_audit/cleanup/api-ownership.csv"),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }
    inventory <- cleanup.namespace.inventory(root)
    errors <- character()

    if (anyDuplicated(ledger$key)) {
        errors <- c(errors, "Ownership ledger contains duplicate keys.")
    }
    missing <- setdiff(inventory$key, ledger$key)
    extra <- setdiff(ledger$key, inventory$key)
    if (length(missing)) {
        errors <- c(errors, paste(
            "Namespace items missing from ledger:",
            paste(missing, collapse = ", ")
        ))
    }
    if (length(extra)) {
        errors <- c(errors, paste(
            "Ledger items absent from namespace:",
            paste(extra, collapse = ", ")
        ))
    }

    required <- c(
        "key", "item_type", "item_name", "source_file", "ownership",
        "intended_package", "current_status", "direct_callers",
        "transitive_callers", "tests", "documentation",
        "known_downstream_callers", "lifecycle_action", "target_release",
        "protected_call_graph"
    )
    absent.columns <- setdiff(required, names(ledger))
    if (length(absent.columns)) {
        errors <- c(errors, paste(
            "Ledger columns missing:", paste(absent.columns, collapse = ", ")
        ))
    } else {
        for (column in required) {
            if (any(is.na(ledger[[column]]) | !nzchar(ledger[[column]]))) {
                errors <- c(errors, paste("Blank ledger field:", column))
            }
        }
    }

    allowed <- c(
        "CORE-ANALYSIS", "CORE-ASSOCIATION", "CORE-PRIVATE", "DGRAPHS",
        "GFLOWX", "DOMAIN", "EXAMPLE", "DELETE", "PROTECTED"
    )
    bad.labels <- setdiff(unique(ledger$ownership), allowed)
    if (length(bad.labels)) {
        errors <- c(errors, paste(
            "Unknown ownership labels:", paste(bad.labels, collapse = ", ")
        ))
    }

    definitions <- cleanup.function.definitions(root)
    graph <- cleanup.call.graph(definitions)
    watched <- ledger$ownership %in% c("DGRAPHS", "GFLOWX", "DOMAIN")
    for (i in which(watched)) {
        actual <- cleanup.collapse(graph$reverse[[ledger$item_name[[i]]]])
        baseline <- ledger$direct_callers[[i]]
        baseline.set <- if (baseline == "none") character() else
            strsplit(baseline, ";", fixed = TRUE)[[1L]]
        actual.set <- if (actual == "none") character() else
            strsplit(actual, ";", fixed = TRUE)[[1L]]
        new.callers <- setdiff(actual.set, baseline.set)
        if (length(new.callers)) {
            errors <- c(errors, paste0(
                ledger$key[[i]], " gained internal callers: ",
                paste(new.callers, collapse = ", ")
            ))
        }
    }
    errors
}

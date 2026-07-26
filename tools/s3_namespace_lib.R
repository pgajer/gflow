s3.split.values <- function(x) {
    if (is.na(x) || !nzchar(x)) {
        return(character())
    }
    trimws(strsplit(x, ";", fixed = TRUE)[[1L]])
}

s3.read.ownership <- function(root) {
    read.csv(
        file.path(root, "split_audit/cleanup/s3-class-ownership.csv"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

s3.validate.namespace <- function(root) {
    ownership <- s3.read.ownership(root)
    inventory <- cleanup.namespace.inventory(root)
    ledger <- read.csv(
        file.path(root, "split_audit/cleanup/api-ownership.csv"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    definitions <- cleanup.function.definitions(root)
    graph <- cleanup.call.graph(definitions)
    errors <- character()

    required <- c(
        "class", "ownership", "constructor", "validator", "stable_fields",
        "methods", "documentation", "tests", "status"
    )
    missing.columns <- setdiff(required, names(ownership))
    if (length(missing.columns)) {
        return(paste(
            "S3 ownership columns missing:",
            paste(missing.columns, collapse = ", ")
        ))
    }

    if (anyDuplicated(ownership$class)) {
        errors <- c(errors, "S3 ownership classes must be unique.")
    }
    if (any(ownership$status != "retained")) {
        errors <- c(errors, "S3 ownership contains a non-retained class.")
    }

    registered <- inventory[inventory$item_type == "s3", , drop = FALSE]
    protected.keys <- ledger$key[
        ledger$item_type == "s3" & ledger$ownership == "PROTECTED"
    ]
    in.scope <- registered[!registered$key %in% protected.keys, , drop = FALSE]

    expected.methods <- unlist(
        lapply(ownership$methods, s3.split.values),
        use.names = FALSE
    )
    missing.methods <- setdiff(in.scope$item_name, expected.methods)
    stale.methods <- setdiff(expected.methods, in.scope$item_name)
    if (length(missing.methods)) {
        errors <- c(errors, paste(
            "Registered in-scope S3 methods lack class ownership:",
            paste(missing.methods, collapse = ", ")
        ))
    }
    if (length(stale.methods)) {
        errors <- c(errors, paste(
            "S3 ownership lists unregistered methods:",
            paste(stale.methods, collapse = ", ")
        ))
    }

    in.scope.ledger <- ledger[
        ledger$key %in% in.scope$key,
        ,
        drop = FALSE
    ]
    undocumented <- in.scope.ledger$item_name[
        in.scope.ledger$documentation == "none"
    ]
    if (length(undocumented)) {
        errors <- c(errors, paste(
            "Registered in-scope S3 methods lack documentation:",
            paste(undocumented, collapse = ", ")
        ))
    }

    for (i in seq_len(nrow(ownership))) {
        row <- ownership[i, ]
        methods <- s3.split.values(row$methods)
        docs <- s3.split.values(row$documentation)
        tests <- s3.split.values(row$tests)

        absent.symbols <- setdiff(
            c(row$constructor, row$validator, methods),
            names(definitions)
        )
        if (length(absent.symbols)) {
            errors <- c(errors, paste0(
                row$class, " references undefined symbols: ",
                paste(absent.symbols, collapse = ", ")
            ))
        }

        missing.docs <- docs[!file.exists(file.path(root, docs))]
        if (length(missing.docs)) {
            errors <- c(errors, paste0(
                row$class, " references missing documentation: ",
                paste(missing.docs, collapse = ", ")
            ))
        }
        missing.tests <- tests[!file.exists(file.path(root, tests))]
        if (length(missing.tests)) {
            errors <- c(errors, paste0(
                row$class, " references missing tests: ",
                paste(missing.tests, collapse = ", ")
            ))
        }

        for (method in intersect(methods, names(graph$forward))) {
            if (!row$validator %in% graph$forward[[method]]) {
                errors <- c(errors, paste0(
                    method, " does not invoke ", row$validator, "."
                ))
            }
        }
    }

    unique(errors)
}

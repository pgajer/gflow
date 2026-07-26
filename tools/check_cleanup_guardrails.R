#!/usr/bin/env Rscript

script.arg <- grep("^--file=", commandArgs(), value = TRUE)
script.path <- if (length(script.arg)) {
    normalizePath(sub("^--file=", "", script.arg[[1L]]))
} else {
    normalizePath("tools/check_cleanup_guardrails.R")
}
root <- dirname(dirname(script.path))
source(file.path(root, "tools/cleanup_ledger_lib.R"))

errors <- c(
    cleanup.validate.ledger(root),
    cleanup.check.protected.surface(root)
)

retired.path <- file.path(
    root, "split_audit/cleanup/retired-exports.txt"
)
retired <- readLines(retired.path, warn = FALSE)
retired <- retired[nzchar(trimws(retired)) & !grepl("^#", trimws(retired))]
current <- cleanup.namespace.inventory(root)
reappeared <- intersect(retired, current$item_name[current$item_type == "export"])
if (length(reappeared)) {
    errors <- c(errors, paste(
        "Retired exports reappeared:", paste(reappeared, collapse = ", ")
    ))
}

if (length(errors)) {
    cat(paste0("ERROR: ", errors, collapse = "\n"), "\n")
    quit(status = 1L)
}

cat(
    "Cleanup guardrails pass for ",
    sum(current$item_type == "export"), " exports and ",
    sum(current$item_type == "s3"), " S3 methods.\n",
    sep = ""
)

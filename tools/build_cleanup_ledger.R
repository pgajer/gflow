#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
scan.downstream <- !("--no-downstream" %in% args)
script.arg <- grep("^--file=", commandArgs(), value = TRUE)
script.path <- if (length(script.arg)) {
    normalizePath(sub("^--file=", "", script.arg[[1L]]))
} else {
    normalizePath("tools/build_cleanup_ledger.R")
}
root <- dirname(dirname(script.path))

source(file.path(root, "tools/cleanup_ledger_lib.R"))

ledger <- cleanup.build.ledger(root, scan.downstream = scan.downstream)
ledger.path <- file.path(
    root, "split_audit/cleanup/api-ownership.csv"
)
write.csv(ledger, ledger.path, row.names = FALSE, na = "")

surface <- cleanup.build.protected.surface(root)
surface.path <- file.path(
    root, "split_audit/cleanup/protected-basin-surface.txt"
)
writeLines(surface, surface.path, useBytes = TRUE)

cat(
    "Wrote", nrow(ledger), "ownership rows (",
    sum(ledger$item_type == "export"), " exports, ",
    sum(ledger$item_type == "s3"), " S3 methods).\n",
    "Protected surface: ",
    sum(grepl("^FILE\\t", surface)), " files, ",
    sum(grepl("^SYMBOL\\t", surface)), " mixed-file symbols, ",
    sum(grepl("^NATIVE\\t", surface)), " native registrations.\n",
    sep = ""
)

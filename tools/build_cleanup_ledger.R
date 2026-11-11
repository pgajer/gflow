#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
scan.downstream <- "--scan-downstream" %in% args
script.arg <- grep("^--file=", commandArgs(), value = TRUE)
script.path <- if (length(script.arg)) {
    normalizePath(sub("^--file=", "", script.arg[[1L]]))
} else {
    normalizePath("tools/build_cleanup_ledger.R")
}
root <- dirname(dirname(script.path))
output.root <- Sys.getenv("GFLOW_LEDGER_OUTPUT", unset = root)
output.root <- normalizePath(output.root, mustWork = TRUE)
if (scan.downstream && (identical(output.root, root) ||
    startsWith(output.root, paste0(root, "/")))) {
    stop("For downstream scans, set GFLOW_LEDGER_OUTPUT to a directory outside the source checkout.")
}

source(file.path(root, "tools/cleanup_ledger_lib.R"))

ledger <- cleanup.build.ledger(root, scan.downstream = scan.downstream)
output.dir <- file.path(output.root, "split_audit/cleanup")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
ledger.path <- file.path(output.dir, "api-ownership.csv")
write.csv(ledger, ledger.path, row.names = FALSE, na = "")

surface <- cleanup.build.protected.surface(root)
surface.path <- file.path(
    output.dir, "protected-basin-surface.txt"
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

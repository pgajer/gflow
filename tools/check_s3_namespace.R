args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
script.path <- if (length(file.arg)) {
    sub("^--file=", "", file.arg[[1L]])
} else {
    normalizePath("tools/check_s3_namespace.R")
}
root <- normalizePath(file.path(dirname(script.path), ".."))

source(file.path(root, "tools/cleanup_ledger_lib.R"))
source(file.path(root, "tools/s3_namespace_lib.R"))

errors <- s3.validate.namespace(root)
if (length(errors)) {
    cat(paste0("ERROR: ", errors, "\n"), sep = "")
    quit(status = 1L)
}

ownership <- s3.read.ownership(root)
cat(
    "S3 namespace ownership passes for",
    nrow(ownership), "in-scope classes and",
    sum(lengths(lapply(ownership$methods, s3.split.values))),
    "methods.\n"
)

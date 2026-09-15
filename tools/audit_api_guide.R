#!/usr/bin/env Rscript
# Run from the package root. Optional --installed also checks the loaded package.
ns <- parse("NAMESPACE")
calls <- function(name) Filter(function(x) is.call(x) &&
    identical(x[[1L]], as.name(name)), as.list(ns))
exports <- vapply(calls("export"), function(x) as.character(x[[2L]]), "")
imports <- calls("importFrom")
import.names <- vapply(imports, function(x) as.character(x[[3L]]), "")
reexports <- intersect(exports, import.names)
s3 <- vapply(calls("S3method"), function(x)
    paste(as.character(x[[2L]]), as.character(x[[3L]]), sep = ","), "")
guide <- readLines("vignettes/function-guide.Rmd", warn = FALSE)
section <- function(name) {
    start <- which(guide == paste0("<!-- ", name, "-start -->"))
    end <- which(guide == paste0("<!-- ", name, "-end -->"))
    stopifnot(length(start) == 1L, length(end) == 1L, start < end)
    guide[seq.int(start + 1L, end - 1L)]
}
rows <- grep("^\\| `[^`]+[(][)]` \\|", section("export-catalog"), value = TRUE)
catalog <- sub("^\\| `([^`]+)[(][)]`.*$", "\\1", rows)
check.set <- function(actual, expected, label) {
    missing <- setdiff(expected, actual)
    extra <- setdiff(actual, expected)
    duplicates <- unique(actual[duplicated(actual)])
    if (length(c(missing, extra, duplicates))) stop(sprintf(
        "%s: missing [%s]; extra [%s]; duplicated [%s]", label,
        paste(missing, collapse = ", "), paste(extra, collapse = ", "),
        paste(duplicates, collapse = ", ")))
}
check.set(catalog, exports, "Export catalog")
stopifnot(all(grepl("[|] (Start|Advanced|Archive|Retired|Re-export) [|]$", rows)))
check.set(catalog[grepl("[|] Re-export [|]$", rows)], reexports, "Re-exports")
method.rows <- grep("^\\| `", section("s3-catalog"), value = TRUE)
method.catalog <- unlist(lapply(method.rows, function(row) {
    fields <- regmatches(row, gregexpr("`[^`]+`", row))[[1L]]
    fields <- substring(fields, 2L, nchar(fields) - 1L)
    paste(fields[[1L]], fields[-1L], sep = ",")
}), use.names = FALSE)
check.set(method.catalog, s3, "S3 catalog")

# Parse source without evaluating it, including unevaluated defaults.
definitions <- list()
for (path in list.files("R", "[.]R$", full.names = TRUE)) {
    for (expr in as.list(parse(path))) {
        if (!is.call(expr) || !as.character(expr[[1L]])[1L] %in% c("<-", "=") ||
            length(expr) != 3L || !is.call(expr[[3L]]) ||
            !identical(expr[[3L]][[1L]], as.name("function"))) next
        definitions[[as.character(expr[[2L]])]] <- expr[[3L]]
    }
}
stopifnot(all(setdiff(exports, reexports) %in% names(definitions)))
method.names <- sub(",", ".", s3, fixed = TRUE)
stopifnot(all(method.names %in% names(definitions)))
retired <- catalog[grepl("[|] Retired [|]$", rows)]
stopifnot(length(retired) == 5L, all(vapply(retired, function(name)
    ".stop.retired.basin.function" %in% all.names(definitions[[name]]), logical(1))))

# Check installed-help coverage against the generated Rd aliases.
aliases <- unlist(lapply(list.files("man", "[.]Rd$", full.names = TRUE), function(p) {
    rd <- tools::parse_Rd(p)
    vapply(Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"), rd),
           function(x) paste(as.character(x), collapse = ""), "")
}), use.names = FALSE)
stopifnot(all(exports %in% aliases), "gflow-migration" %in% aliases)

if ("--installed" %in% commandArgs(TRUE)) {
    loadNamespace("gflow")
    check.set(getNamespaceExports("gflow"), exports, "Installed exports")
    stopifnot(identical(gflow::vertices, dgraphs::vertices))
    for (method in s3) {
        parts <- strsplit(method, ",", fixed = TRUE)[[1L]]
        stopifnot(is.function(utils::getS3method(parts[1L], parts[2L],
            optional = TRUE, envir = asNamespace("gflow"))))
    }
    for (name in retired) {
        error <- tryCatch(do.call(getExportedValue("gflow", name), list()),
                          error = identity)
        stopifnot(inherits(error, "gflow_basin_lifecycle_error"))
    }
    v <- utils::vignette(package = "gflow")$results
    stopifnot(all(c("function-guide", "example-graphs-and-fields") %in% v[, "Item"]))
}
cat(sprintf(paste0("Guide coverage: %d exports (%d local functions, %d re-export); ",
                   "%d S3 registrations; %d exported method names; %d retirement stubs.\n"),
            length(exports), length(setdiff(exports, reexports)), length(reexports),
            length(s3), length(intersect(exports, method.names)), length(retired)))
cat("No missing or duplicate catalog rows; generated help aliases verified.\n")

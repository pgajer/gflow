args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
script.path <- if (length(file.arg)) {
    sub("^--file=", "", file.arg[[1L]])
} else {
    normalizePath("tools/check_phase7_ownership.R")
}
root <- normalizePath(file.path(dirname(script.path), ".."))

dependency.path <- file.path(
    root, "split_audit/cleanup/dependency-ownership.csv"
)
native.path <- file.path(
    root, "split_audit/cleanup/native-symbol-ownership.csv"
)
if (!file.exists(dependency.path) || !file.exists(native.path)) {
    stop("Phase 7 ownership ledgers are missing; run tools/build_phase7_ownership.R.")
}

dependency <- read.csv(
    dependency.path,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
native <- read.csv(
    native.path,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

description <- read.dcf(file.path(root, "DESCRIPTION"))[1L, ]
split.dependencies <- function(value) {
    if (is.na(value) || !nzchar(value)) return(character())
    trimws(sub("\\s*\\(.*$", "", strsplit(value, ",", fixed = TRUE)[[1L]]))
}
declared.keys <- unlist(lapply(
    c("Depends", "Imports", "Suggests", "LinkingTo"),
    function(field) {
        paste(field, split.dependencies(description[[field]]), sep = ":")
    }
))
declared.keys <- c(
    declared.keys,
    paste(
        "VignetteBuilder",
        split.dependencies(description[["VignetteBuilder"]]),
        sep = ":"
    ),
    paste(
        "SystemRequirements",
        c("C++17", "GNU make", "OpenMP"),
        sep = ":"
    )
)
declared.keys <- declared.keys[nzchar(sub("^[^:]+:", "", declared.keys))]

errors <- character()
if (!setequal(dependency$key, declared.keys)) {
    errors <- c(
        errors,
        paste(
            "Dependency ownership coverage differs from DESCRIPTION:",
            paste(
                c(
                    paste0("missing=", setdiff(declared.keys, dependency$key)),
                    paste0("stale=", setdiff(dependency$key, declared.keys))
                ),
                collapse = ", "
            )
        )
    )
}
if (any(!nzchar(dependency$ownership)) ||
    any(!nzchar(dependency$feature)) ||
    any(!nzchar(dependency$policy))) {
    errors <- c(errors, "Dependency ownership contains blank required fields.")
}

init.lines <- readLines(file.path(root, "src/init.c"))
registered <- sub(
    "^[[:space:]]*\\{\"([^\"]+)\".*",
    "\\1",
    grep("^[[:space:]]*\\{\"[^\"]+\"", init.lines, value = TRUE)
)
if (anyDuplicated(registered)) {
    errors <- c(errors, "src/init.c contains duplicate native registrations.")
}
if (!setequal(native$symbol, registered)) {
    errors <- c(
        errors,
        paste(
            "Native ownership coverage differs from src/init.c:",
            paste(
                c(
                    paste0("missing=", setdiff(registered, native$symbol)),
                    paste0("stale=", setdiff(native$symbol, registered))
                ),
                collapse = ", "
            )
        )
    )
}
if (any(native$r_callers == "none")) {
    errors <- c(
        errors,
        paste(
            "Registered native symbols without a live R/test caller:",
            paste(native$symbol[native$r_callers == "none"], collapse = ", ")
        )
    )
}
if (any(native$native_sources == "none")) {
    errors <- c(
        errors,
        paste(
            "Registered native symbols without a source owner:",
            paste(native$symbol[native$native_sources == "none"], collapse = ", ")
        )
    )
}

if (length(errors)) {
    cat(paste0("ERROR: ", errors, "\n"), sep = "")
    quit(status = 1L)
}

cat(
    "Phase 7 ownership passes for", nrow(dependency),
    "dependency declarations and", nrow(native),
    "native registrations.\n"
)

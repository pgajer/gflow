active.major.minor <- paste(
    R.version$major,
    strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][[1L]],
    sep = "."
)

library.paths <- normalizePath(
    .libPaths(),
    winslash = "/",
    mustWork = TRUE
)
version.pattern <- "(?<=/)[0-9]+\\.[0-9]+(?=[/-])"
library.versions <- unlist(
    regmatches(
        library.paths,
        gregexpr(version.pattern, library.paths, perl = TRUE)
    ),
    use.names = FALSE
)
foreign.versions <- setdiff(unique(library.versions), active.major.minor)
if (length(foreign.versions) > 0L) {
    stop(
        "Library paths from another R minor version are active: ",
        paste(foreign.versions, collapse = ", "),
        ". Active R minor version: ",
        active.major.minor,
        call. = FALSE
    )
}

if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp is not installed for the active R toolchain.", call. = FALSE)
}

rcpp.description <- packageDescription("Rcpp")
rcpp.built.version <- sub(
    "^R ([0-9]+\\.[0-9]+).*$",
    "\\1",
    rcpp.description$Built
)
if (!identical(rcpp.built.version, active.major.minor)) {
    stop(
        "Rcpp was built for R ",
        rcpp.built.version,
        " but the active R minor version is ",
        active.major.minor,
        ".",
        call. = FALSE
    )
}

if (!dir.exists(R.home("include"))) {
    stop("The active R header directory does not exist.", call. = FALSE)
}

compiled.value <- Rcpp::evalCpp("1 + 1")
if (!identical(compiled.value, 2L)) {
    stop("The Rcpp compilation probe returned an unexpected value.", call. = FALSE)
}

cat(
    "R toolchain preflight passed\n",
    "  R: ", R.version.string, "\n",
    "  R home: ", R.home(), "\n",
    "  R headers: ", R.home("include"), "\n",
    "  Rcpp: ", rcpp.description$Version, " (built for R ",
    rcpp.built.version, ")\n",
    "  Rcpp library: ", find.package("Rcpp"), "\n",
    sep = ""
)
